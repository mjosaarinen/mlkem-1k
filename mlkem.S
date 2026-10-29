//  mlkem.S
//  === A size-oriented ML-KEM (FIPS 203) for x86-64 (SysV ABI).
//
//  No precomputed tables or matrix storage; encoding and Keccak are shared
//  across all three FIPS 203 parameter sets.
//
//  NOT for production use: variable time (DIV, branches on secret data),
//  no wiping of secrets, no input validation (FIPS 203 Sect. 7.2/7.3
//  modulus / hash checks are left to the caller).
//
//  Build (GCC preprocesses .S files):
//      gcc -c -DMLKEM_K=2 mlkem.S          // ML-KEM-512  (default)
//      gcc -c -DMLKEM_K=3 mlkem.S          // ML-KEM-768
//      gcc -c -DMLKEM_K=4 mlkem.S          // ML-KEM-1024
//
//  C interface (the deterministic "internal" functions of FIPS 203):
//
//  void mlkem_keygen(uint8_t ek[EK], uint8_t dk[DK],
//                    const uint8_t d[32], const uint8_t z[32]);   // Alg. 16
//  void mlkem_encaps(uint8_t ct[CT], uint8_t k[32],
//                    const uint8_t ek[EK], const uint8_t m[32]);  // Alg. 17
//  void mlkem_decaps(uint8_t k[32], const uint8_t ct[CT],
//                    const uint8_t dk[DK]);                       // Alg. 18
//
//  Outputs must not overlap inputs. No static data is used: everything
//  lives on the stack (peak use is about 15 kB for ML-KEM-1024 decaps).
//
//  Structure: the exported functions share one prologue (frame_enter).
//  K-PKE.Encrypt is inlined into encaps and K-PKE.Decrypt into decaps;
//  decaps re-encrypts by calling encaps, and both share the row product
//  "dot". All polynomial products are computed in the normal domain by a
//  negacyclic schoolbook multiplication; the O(n^2) NTT evaluator only
//  converts the NTT-domain values that FIPS 203 defines (A^, t^, s^).
//  SampleNTT, SamplePolyCBD and ByteDecode_12 of t^ / s^ share a sampler.
//
//  Internal routines use custom register contracts, documented at each
//  routine; callers rely on them exactly (e.g. returned start / end
//  pointers in RCX / RDI). Polynomials are 256 x uint16 (512 bytes).
//  Size tricks used throughout: push imm8 / pop reg for small constants
//  and push reg / pop reg for register copies, string instructions,
//  fall-through between routines, and values that routines are known to
//  leave behind (e.g. RCX = 0 or EDX bits 16+ after Keccak, SF and ZF
//  selecting the sampler mode), each noted where it is relied on.

    .intel_syntax noprefix

//  --- parameter sets (FIPS 203, Table 2)

#ifndef MLKEM_K
#define MLKEM_K     2
#endif

#if     MLKEM_K == 2
#define ETA1        3
#define ETA2        2
#define DU          10
#define DV          4
#elif   MLKEM_K == 3
#define ETA1        2
#define ETA2        2
#define DU          10
#define DV          4
#elif   MLKEM_K == 4
#define ETA1        2
#define ETA2        2
#define DU          11
#define DV          5
#else
#error  "MLKEM_K must be 2, 3 or 4."
#endif

#define Q           3329
#define POLY        512                         // bytes per polynomial
#define EK_BYTES    (384*MLKEM_K+32)
#define DK_BYTES    (768*MLKEM_K+96)
#define CT_BYTES    (32*(DU*MLKEM_K+DV))

//  Keccak sponge "modes": rate in bytes | (domain padding byte << 8).
#define SHA3_256    0x0688                      // H: rate 136, pad 0x06
#define SHA3_512    0x0648                      // G: rate  72, pad 0x06
#define SHAKE128    0x1FA8                      // XOF: rate 168, pad 0x1F
#define SHAKE256    0x1F88                      // PRF, J: rate 136, pad 0x1F

//  Frame layout shared by the three bodies and dot. The frame starts with a
//  HDR-byte header (m || H(ek) || K || r in encaps, d || k || rho || sigma
//  in keygen, z || c in decaps) at RBP-HDR; RBP points just after it, at
//  the vector (s, y or u), then come the accumulators (e1[K] -> u[K] and
//  e2 -> v in encaps, v -> w in decaps, e[K] -> t[K] in keygen), then
//  dot's scratch (R9). HDR = 128 keeps all polynomials at 0 mod 8
//  (ntt_common needs bit 1 clear), the header within disp8 reach of RBP
//  and fits m || H(ek) || K || r, so G writes K || r right after its
//  input.
#define HDR         128
#define FUNC(f)     .type f, @function; f:
#define ENDF(f)     .size f, . - f

    .text

//  ========================================================================
//  Arithmetic mod q (the evaluator below the sampler, see ntt_common)
//  ========================================================================

//  The NTT is computed directly from its definition as an O(n^2)
//  evaluation (FIPS 203 eq. 4.8 / Alg. 9-10 give the same values):
//      fwd:  NTT[2i+b]  = sum_j f[2j+b] * z^(o(i)*j)
//      inv:  f[2j+b]    = 128^-1 * sum_i NTT[2i+b] * z^(-o(i)*j)
//  with z = 17 and o(i) = 2*BitRev7(i)+1. Exponents are taken mod 256 and
//  looked up in a table T[e] = c*g^e (fwd: c = 1, g = 17; inv:
//  c = 128^-1 = 3303, or -128^-1 for decaps' -s, g = 17^-1) that is
//  pushed onto the stack as P_k = c*m^(k+1), m = g^-1, so T[e] = P_255-e
//  sits at [rsp+8e].
//  polymul_acc shares the output path: every output is pushed as an
//  unreduced qword and the final loop pops, reduces and stores them to
//  RSI[255..0]. All entry points need R8 = q (sample_cbd's exit sets it
//  before invntt, which precedes every polymul_acc and ntt in dot and
//  keygen), preserve RSI, RBX, RBP, R8-R15 and return fully reduced
//  coefficients in [0, q). They clobber RAX, RCX, RDX, RDI (transforms
//  keep it) and use about 4 kB of stack.

//  poly_decode: Decompress_d(ByteDecode_d) of the ciphertext (decaps):
//  BH polynomials with d = ECX (exactly; DU), then one with d = DV (BH = 0:
//  just that one), stored consecutively from RDI. The input is read from
//  RSI (kept) at the running bit position EDX; values are below q. BL = 0
//  counts the 256 coefficients of each polynomial. Then falls through into
//  dot with RDI = R9 = end of the last polynomial (dot's accumulator ends
//  there, its scratch starts there), BL = K and BH = 0xFF.

FUNC(poly_decode)
    xor     edx, edx
.Ldec_coef:
    xor     eax, eax                    // get_bits keeps the bits above d
    call    get_bits
    imul    eax, eax, Q
    shr     eax, cl                     // carry = original product bit d-1
    adc     eax, 0                      // Decompress_d rounding
    stosw
    inc     bl                          // 256 coefficients from BL = 0
    jnz     .Ldec_coef
    dec     bh
    jg      .Ldec_coef                  // BH > 0: same d
    mov     cl, DV
    jz      .Ldec_coef                  // BH = 0: v
    mov     bl, MLKEM_K                 // dot's decode row: -s in decaps
    mov     r9, rdi                     // dot's scratch
    //  Fall through into dot.
ENDF(poly_decode)

//  dot: the row product shared by K-PKE.KeyGen, K-PKE.Encrypt and
//  K-PKE.Decrypt, in the normal domain. RDI = end of the accumulator
//  (all callers have just written it or got it back from dot), vec[j] at
//  RBP+512j, R9 = scratch polynomial (kept: cbd_vec and poly_decode set it
//  to the free space after their last polynomial, which is the end of the
//  first accumulator), R10 = rho, ek or dk.
//  acc += sum_j M[j] * vec[j], j = 0..K-1 (BH), with
//      BL < K:  M[j] = NTT^-1(SampleNTT(rho || BL || j)) = A[j][BL],
//               rho = R10 (encaps rows)
//      BL = ~i: M[j] = NTT^-1(SampleNTT(rho || j || i)) = A[i][j] (keygen)
//      BL = K:  M[j] = +-NTT^-1(ByteDecode_12(R10 + 384j)), t (EBX bits
//               16+ = 0, encaps) or -s (EBX bits 16+ = DEC_NEG, decaps;
//               see sample_cbd and invntt); R10 then ends 384K bytes
//               further (rho, or ek inside dk).
//  All rows run through sample_cbd with the flags of CMP AL, K (AL = BL;
//  the keygen swap keeps them): ZF = 1 iff BL = K, otherwise SF = 1 (BL - K
//  wraps negative for BL < K and for BL = ~i).
//  dot keeps EBX bits 16+ (keygen and encaps keep them 0, as frame_enter
//  sets them). Needs RCX < 256 (sample_cbd; callers leave 0 or DV).
//  Returns BH = K, RDI = RSI = accumulator and RCX = 0.

FUNC(dot)
    lea     rax, [rdi-POLY]
    push    rax                         // accumulator
    mov     bh, 0
.Ldot_term:
    push    r10
    pop     rsi                         // rho, or the next 384-byte block
    mov     eax, ebx
    cmp     al, MLKEM_K                 // ZF = 1: decode
    jbe     .Ldot_seed                  // encaps: BL = i < K, seed i || j
    xchg    al, ah                      // keygen: BL = ~i > K, seed j || i
    not     ah
.Ldot_seed:
    call    sample_cbd                  // R9 = A[j][i], t[j] or -s[j]
    pop     rsi
    push    rsi                         // accumulator
    push    rbp
    pop     rcx                         // a = vec[BH] (polymul_acc reads BH)
    call    polymul_acc
    xor     ecx, ecx
    inc     bh
    cmp     bh, MLKEM_K
    jb      .Ldot_term
    pop     rdi                         // RDI = RSI = accumulator
.Ldot_ret:
    ret
ENDF(dot)

//  ========================================================================
//  Sampling
//  ========================================================================


//  sample_cbd is the sampler body shared by three modes; the flags at
//  entry select one (no instruction before the JS / JE changes them):
//  SF = 1 SampleNTT, ZF = 1 ByteDecode_12 of a key polynomial, SF = ZF = 0
//  SamplePolyCBD. All modes need RCX < 256 and start with EBX = 0.
//
//  The output polynomial is at R9 (kept) in all modes.
//
//  SampleNTT(rho || j || i) (Algorithm 7): RSI = rho (32 bytes), AL = j,
//  AH = i.
//  SHAKE128 is squeezed one block at a time; 12-bit candidates are read
//  from the state (1344 = 112 * 12 bits per block, so they never straddle
//  blocks) and accepted (get_bits' CF = 1) below the threshold R8W = q.
//
//  Decode: R10 = RSI = a 384-byte ByteEncode_12 block (t^[j] in ek or
//  s^[j] in dk; the seed copy below reads its first 32 bytes, unused).
//  The same candidate loop reads 256 12-bit values from R10 with
//  R8W = SHAKE256 = 0x1F88 > 4095 as the threshold, so every value is
//  accepted: 3329..4095 stay and invntt reduces them mod q (FIPS 203 uses
//  ByteDecode_12's output mod q, unlike SampleNTT's candidates). R10
//  advances by 384.
//
//  Both modes fall through into invntt, so they return NTT^-1 of the
//  sampled / decoded polynomial in the normal domain at R9 (and RSI). R8
//  holds q in its upper half (set on entry) and the threshold in R8W; the
//  exit shifts q down for ntt_common (SampleNTT's R8W is q as well).
//  get_bits keeps EAX bits 12+, which give invntt's c (see there): in
//  SampleNTT they are 0 (keccak_f1600 returns EAX = 127, later EAX is the
//  previous candidate); in decode they are dot's EBX bits 12+ (EAX = EBX
//  at entry; BH = j < 16, so bits 12..15 are 0).
//
//  cbd_vec: SamplePolyCBD_eta(PRF_eta(s, N)) (Algorithm 8 and Sect. 4.1)
//  for N = BL..BH-1, stored consecutively from RBP (the frame vector).
//  RSI = s (kept); BL starts at 0 (still frame_enter's EBX = 0 in both
//  callers), RCX < 256.
//  eta = eta1, except eta2 for N >= K when BH = 2K+1 (encaps' e1 and e2;
//  keygen has BH = 2K): only ML-KEM-512 has eta1 != eta2, and sample_cbd
//  picks eta from AX = BH:N. Ends with BL = BH, RCX = 0 and RDI = R9 after
//  the last polynomial (dot's scratch). Each pass sets R9 = RDI and
//  enters sample_cbd with the flags of CMP BH, BL: BH > BL gives
//  SF = ZF = 0 (CBD); the CBD exit loops back. Each coefficient runs two
//  halves of "negate, then add 1 - b for each of the next eta bits b"
//  from q: -(-q + eta - x) + eta - y = q + x - y, in [q - eta, q + eta]
//  (not reduced; polymul_acc accepts that, see there). The value is
//  negative after the first half only, so its sign ends the coefficient.
//
//  PRF_eta needs 64 * eta <= 192 bytes, more than the SHAKE256 rate (136)
//  for eta = 3: the first block is copied to RSP, just below the state,
//  and permuting the state puts the second block right after it.
//
//  Frame: 136 bytes, then the Keccak state at RSP+136. The seed (32 bytes
//  || 2-byte suffix) ends at the state. R9 = polynomial (Keccak keeps it).

FUNC(cbd_vec)
    push    rbp
    pop     rdi                         // the vector, after the header
.Lcv_poly:
    mov     r9, rdi                     // polynomial; at the end dot's scratch
    cmp     bh, bl                      // BH > BL: SF = ZF = 0 selects CBD
    jbe     .Ldot_ret
    push    rsi                         // s
    mov     eax, ebx                    // N = BL, AH = BH
ENDF(cbd_vec)

FUNC(sample_cbd)
    push    rbx
    push    0
    pop     rbx                         // EBX = 0 (flags kept)
    mov     r8d, SHAKE256 | Q << 16     // q for ntt_common, see .Lsn_done
    enter   336, 0                      // 136 + 200 bytes
    lea     rdi, [rsp+136-34]           // the seed ends at the state
    push    rdi
    mov     cl, 32                      // RCX < 256 at entry
    rep     movsb
    stosw                               // seed = s || N, or rho || j || i
    pop     rsi
    push    34
    pop     rdx
    js      .Lsn_absorb                 // no flags changed since entry
    je      .Lsn_decode

    dec     edx                         // s || N
#if ETA1 != ETA2
    //  eta = eta2 for encaps' e1 / e2 (AX = BH:N = 2K+1:N, N >= K), else
    //  eta1 (keygen has BH = 2K): CF = 1 iff AX >= (2K+1):K.
    add     ax, -((2*MLKEM_K+1) << 8 | MLKEM_K)
    sbb     ebx, -ETA1                  // EBX = eta1 - CF (eta2 = eta1 - 1)
#else
    mov     bl, ETA1                    // eta1 = eta2
#endif
    call    keccak_absorb               // R8W = SHAKE256: rate 136
    push    rdi
    pop     rsi                         // state
    push    rsp
    pop     rdi
    mov     cl, 136                     // RCX = 0 after Keccak
    rep     movsb                       // first block to RSP, RDI = state
    call    keccak_f1600                // second block right after it
    push    r9
    pop     rdi                         // polynomial
    xor     esi, esi                    // bit position
    //  DL = 0 (Keccak's chi leaves it) counts the coefficients. EAX bits
    //  16+ are 0: keccak_f1600 returns EAX = 127, later EAX is the
    //  previous coefficient.
.Lcbd_coef:
    mov     ax, Q                       // EAX = q
.Lcbd_half:
    neg     eax                         // -q, then q - eta + sum(x)
    mov     ecx, ebx
.Lcbd_bit:
    bt      dword ptr [rsp], esi
    inc     esi                         // keeps CF
    sbb     eax, -1                     // adds 1 - bit, SF = sign
    loop    .Lcbd_bit
    js      .Lcbd_half                  // < 0 after the x bits only
    stosw                               // q + sum(x) - sum(y)
    inc     dl                          // wrap after coefficient 255
    jnz     .Lcbd_coef
    leave                               // RCX = 0 (LOOP), RDI = end
    pop     rbx
    pop     rsi                         // s
    inc     ebx
    jmp     .Lcv_poly

.Lsn_absorb:
    mov     r8b, SHAKE128 & 0xFF        // same padding byte as SHAKE256
    call    keccak_absorb               // EBX = 0 counts accepted coefficients
    mov     r8w, Q                      // threshold; R8D = q * 0x10001
    //  Keccak returns RSI = state + 200 and leaves EDX bits 16+ alone; they
    //  are all ones here (EDX = -1 after absorbing, -256 after a block).
.Lsn_block:
    mov     dx, -1600                   // EDX = -1600: bit 0 of the state
.Lsn_twelve:
    mov     cl, 12                      // RCX = 0 after Keccak; get_bits keeps it
.Lsn_candidate:
    call    get_bits                    // CF = 1: accept
    mov     [r9+rbx*2], ax              // a rejected value is overwritten later
    adc     bl, bh                      // BH = 0; CF after the 256th accept
    jc      .Lsn_done
    cmp     dh, 0xFF                    // EDX = -256 after 1344 = 168*8 bits
    jne     .Lsn_candidate
    call    keccak_f1600
    jmp     .Lsn_block
.Lsn_decode:
    mov     esi, 384
    xadd    r10, rsi                    // RSI = block, R10 = next block
    cdq                                 // EDX = 0: EAX = dot's EBX >= 0
    jmp     .Lsn_twelve                 // R8W = 0x1F88: accept all
.Lsn_done:
    leave
    pop     rbx
    shr     r8d, 16                     // R8 = q (SampleNTT and decode)
    //  Fall through to invntt (R9 = polynomial): SampleNTT returns
    //  A = NTT^-1(A^), decode t = NTT^-1(t^) or -s = -NTT^-1(s^).
ENDF(sample_cbd)

//  ========================================================================
//  NTT and polynomial multiplication (ntt_common)
//  ========================================================================

//  invntt: R9 = polynomial, ntt: RSI = polynomial, transformed in place
//  (both return it in RSI, keep EDI and return RCX = 0, which keygen and
//  dot rely on; inside, RCX = 0 selects the transform and EDI = table
//  multiplier m is also the direction flag). invntt's only entry is
//  sample_cbd's fall-through. It sets AX = 3303 and EDI = m, so
//  c = 3303 + 65536 * (EAX >> 16): sample_cbd leaves EAX >> 16 = 0
//  (c = 128^-1) except in decaps' decode rows, where it is DEC_NEG
//  (c = -128^-1; any 32-bit c keeps DIV in range since m = 17). ntt is
//  the forward entry; its only caller (keygen) passes EAX = c = 1 and
//  EDI = m + 645084q, which has bit 31 set (DIV stays in range since
//  EAX < q is multiplied by EDI < 2^32). Inputs may be up to 4095
//  (unreduced ByteDecode_12).
//  polymul_acc: RSI += a * R9 in Z_q[X]/(X^256+1) with a = RCX + 512*BH
//  (dot passes RCX = RBP and BH = j: a = vec[j]); RSI must not alias a
//  or R9; RCX nonzero selects this mode. b = R9 <= q (b is
//  always a transform output); a and the accumulator may be up to q + 3
//  (CBD samples are left as q + v). The (unused) table is
//  still built from EAX and EDI; EAX may be large on entry, so EDI must be
//  small to keep DIV in range: it is 17 after invntt, which always
//  precedes polymul_acc. Clobbers EDI.

FUNC(invntt)
    mov     ax, 3303                    // c = 128^-1
    push    17
    pop     rdi                         // m = 17
    push    r9
    pop     rsi                         // polynomial
ENDF(invntt)

FUNC(ntt)
    xor     ecx, ecx                    // RCX = 0: transform
ENDF(ntt)

FUNC(ntt_common)
FUNC(polymul_acc)
    push    rbx
    push    rbp                         // RBP restores RSP after the table
    push    rsp
    pop     rbp
    xor     ebx, ebx                    // BL = output index x (and counter)
.Lnt_tab:
    mul     edi                         // product < 2^32
    div     r8d
    xchg    eax, edx
    push    rax
    inc     bl
    jnz     .Lnt_tab
    cdq                                 // EAX < q: EDX = 0 (polymul's i)
    mov     dh, [rbp+9]                 // caller's BH (polymul's j)
    //  x outputs are already pushed, so T[e] is at [rsp+8(e+x)].
.Lnt_out:
    jrcxz   .Lnt_body
    //  polymul: slot = acc[x] + sum_i a[i] * b[x-i] in Z_q[X]/(X^256+1);
    //  terms that wrap around use q - b[x-i]. With acc, a <= q+3 and
    //  b <= q the sum is below q+4 + 256*(q+3)*q < 2^32, so a 32-bit ADD
    //  suffices.
    //  EDX = 256j + i with i = 0 here: the inner loop wraps DL back to 0,
    //  so [RCX+RDX*2] = a[i].
    movzx   eax, word ptr [rsi+rbx*2]
    push    rax
.Lsb_term:
    mov     eax, ebx
    sub     al, dl                      // x - i mod 256, CF = wrapped
    movzx   eax, word ptr [r9+rax*2]
    jnc     1f
    neg     eax
    add     eax, r8d                    // q - b
1:  movzx   edi, word ptr [rcx+rdx*2]
    imul    eax, edi
    add     [rsp], eax
    inc     dl
    jnz     .Lsb_term
.Lnt_next:
    inc     bl
    jnz     .Lnt_out
.Lnt_pop:
    pop     rax
    xor     edx, edx
    div     r8
    dec     bl
    mov     [rsi+rbx*2], dx
    jnz     .Lnt_pop
    leave
    pop     rbx
    ret

    //  Transform output x = 2i+b: slot = sum over n = 2y of in[n+b] * T[e].
    //  RSI is offset by 2b: it toggles 256 times and ends unchanged. This
    //  needs bit 1 of RSI clear; all polynomials are at 0 mod 8 in
    //  frames whose base is 8-aligned (RSP is always a multiple of 8).
    //  fwd: e = o(i)*y; inv (m = 17): e = o(y)*i, so the roles swap.
.Lnt_body:
    push    rcx                         // zero accumulator slot
.Lnt_in:
    mov     eax, ebx
    mov     edx, ecx
    test    edi, edi                    // fwd: EDI >= 2^31
    js      1f
    xchg    eax, edx
1:  mov     ah, 1                       // sentinel: EDX = B'<<9 | ? | o(A)
2:  shr     eax, 1
    rcl     edx, 1
    jnz     2b
    shld    eax, edx, 22                // B = B'/2
    mul     dl
    mov     ah, 0                       // e = o(A)*B mod 256
    add     eax, ebx
    movzx   edx, word ptr [rsi+rcx*2]
    imul    edx, [rsp+rax*8+8]
    add     [rsp], edx                  // 128 terms < 2^31
    add     cl, 2
    jnc     .Lnt_in                     // ends with RCX = 0 again
    xor     rsi, 2
    jmp     .Lnt_next
ENDF(polymul_acc)
ENDF(ntt_common)


//  ========================================================================
//  Encoding, compression (Section 4.2.1)
//  ========================================================================

//  get_bits: replaces the low CL (1..31) bits of EAX by the CL bits at bit
//  position EDX of [RSI], little endian, and keeps EAX's other bits (the
//  callers clear them or rely on them); EDX advances by CL. Ends with
//  CMP AX, R8W: CF = 1 iff AX is below the threshold R8W (sample_cbd).
//  Clobbers nothing else.

FUNC(get_bits)
    push    rcx
.Lgb_bit:
    bt      dword ptr [rsi], edx
    rcr     eax, 1
    inc     edx
    loop    .Lgb_bit
    pop     rcx
    rol     eax, cl                     // top CL collected bits move to the bottom
    cmp     ax, r8w
    ret
ENDF(get_bits)


//  ========================================================================
//  ML-KEM (Section 6), exported
//  ========================================================================

//  frame_enter: shared prologue/epilogue. An exported function starts with
//  "call frame_enter"; the rest of that function runs as a subroutine with
//  RBX and RBP saved, EBX = 0 (bodies change BL / BH, decaps also bits
//  16+), RAX = 1st and R10 = 3rd argument, RSI and RDX still holding the
//  2nd and 3rd, RBP = header + HDR (see the frame layout) for a FRAME-byte
//  scratch area that starts with the header. The 32 bytes at RCX (4th
//  argument: m in encaps, d in keygen, z in decaps, each set up so) are
//  copied to the header first, so bodies start with RDI = header + 32 and
//  RCX = 0. Bodies keep R12-R15.
//  The body's RET returns here, which restores everything (RSP from the
//  slot below the frame) and returns to the caller. The frame base is
//  8-aligned (see ntt_common).

//  FRAME >= HDR + (2K+2)*POLY (vector, accumulators, dot's scratch) and
//  >= 32 + CT_BYTES + 200 (J's input z || c and Keccak state), 8-aligned.
#define FRAME       (1024*MLKEM_K+1152)

FUNC(frame_enter)
    pop     r11                         // body address
    push    rbx
    xor     ebx, ebx                    // bodies get EBX = 0
    xchg    rax, rdi                    // RAX = 1st argument
    mov     r10, rdx                    // caller-saved: R10 is never restored
    enter   FRAME, 0                    // push rbp, RBP = top, allocate
    push    rsp
    pop     rdi                         // RDI = frame base (header)
    push    rbp                         // [base-8] = top
    push    rsi
    push    rcx
    pop     rsi
    push    32
    pop     rcx
    rep     movsb                       // 4th argument into the header
    pop     rsi
    lea     rbp, [rdi+HDR-32]           // RBP = vector, after the header
    call    r11
    pop     rsp                         // back to the top
    pop     rbp
    pop     rbx
    ret
ENDF(frame_enter)


//  mlkem_keygen: ML-KEM.KeyGen_internal(d, z) (Algorithms 16 and 13).
//  RDI = ek, RSI = dk, RDX = d, RCX = z (swapped to RCX = d, RDX = z for
//  frame_enter). Computes t = A*s + e in the normal domain
//  (A = NTT^-1(A^)); s^ = NTT(s) and t^ = NTT(t) are encoded.

//  Header (at RBP-HDR): d || k (33 bytes), then rho || sigma = G(d || k);
//  s[K], then e[K] -> t[K] from RBP.
//  BH = K + n*KG_STEP after n encode steps: below 0x80 for n < K, in
//  [0x80, 0x100) for K <= n < 2K, carry at n = 2K (K = 2, 3, 4).
#define KG_STEP     ((128+MLKEM_K-1)/MLKEM_K)

    .globl  mlkem_keygen
FUNC(mlkem_keygen)
    xchg    rcx, rdx                    // frame_enter copies d
    call    frame_enter
    push    rax                         // ek, popped last
    push    rdx                         // z
    push    rsi                         // dk

    //  (rho, sigma) = G(d || k)
    mov     al, MLKEM_K
    stosb                               // d || k ends at rho
    push    rdi                         // rho
    push    rdi                         // G's output: in place (the state)
    call    hash_G                      // input d || k: 33 bytes

    //  s[i] = CBD(PRF(sigma, i)), e[i] = CBD(PRF(sigma, K+i))
    pop     r10                         // rho, where dot reads it; RSI =
                                        // rho + 32 = sigma
    mov     bh, 2*MLKEM_K               // s[K], e[K] from RBP
    call    cbd_vec

    //  t[i] = e[i] + sum_j A[i][j] * s[j], one dot per row (BL = ~i marks
    //  the transposed seed rho || j || i).
    //  Rows run i = K-1 .. 0, from the last t (cbd_vec ended after it).
    mov     bl, 0x100-MLKEM_K           // BL = ~(K-1); EBX bits 16+ are 0
.Lkg_t:
    call    dot                         // RDI = end of t[i]; dot returns
                                        // RDI = t[i], the end of t[i-1]
    inc     bl
    jnz     .Lkg_t

    //  dk = ByteEncode_12(s^) || ByteEncode_12(t^) || rho || H(ek) || z,
    //  where the middle part ByteEncode_12(t^) || rho is ek.
    pop     rdi                         // dk
    push    rbp
    pop     rsi                         // s, then t
    //  The rows left BH = K and OF = 0 (INC BL from 0xFF). BH steps by
    //  KG_STEP: the K-th step is the only one that crosses 0x7F (OF = 1),
    //  the 2K-th the only one that carries out (see KG_STEP).
.Lkg_enc:
    jno     .Lkg_cursor
    push    rdi                         // ek inside dk, after s^
.Lkg_cursor:
    push    rdi                         // dk cursor (EDI carries m)
    push    1
    pop     rax                         // c = 1
    mov     edi, 1175+645084*Q          // m = 17^-1 (mod q), bit 31 set
    call    ntt                         // s^, t^ (RSI = polynomial)
    pop     rdi
    call    poly_encode                 // RCX = 0 after ntt: d = 12
    add     bh, KG_STEP
    jnc     .Lkg_enc
    push    r10
    pop     rsi                         // rho
    mov     cl, 32                      // poly_encode left RCX = 0
    rep     movsb                       // RDI ends at the H(ek) field
    pop     rsi                         // ek inside dk
    push    rsi
    call    hash_H                      // RDI ends at the z field
    pop     rax                         // ek inside dk
    pop     rsi                         // z
    mov     cl, 32                      // xhash_out left RCX = 0
    rep     movsb
    xchg    rax, rsi
    pop     rdi                         // ek
    mov     cl, EK_BYTES/8              // REP left RCX = 0
    rep     movsq                       // ek = copy from dk
    ret
ENDF(mlkem_keygen)

//  mlkem_encaps: ML-KEM.Encaps_internal(ek, m) (Algorithm 17) with
//  K-PKE.Encrypt (Algorithm 14) inlined, in the normal domain:
//  u[i] = e1[i] + sum_j A[j][i] * y[j] and v = e2 + sum_j t[j] * y[j] + mu,
//  so each accumulator starts as a sampled noise polynomial.
//  RDI = ciphertext, RSI = shared key K, RDX = ek, RCX = m.
//  Returns RDI = end of the ciphertext and RCX = DV (from poly_encode).

//  Frame: header m (0), H(ek) (32), K || r (64, 96) = G(m || H(ek));
//  y[K], e1[K] -> u[K] and e2 -> v from RBP (the order of c = c1 || c2),
//  then the scratch polynomial.

    .globl  mlkem_encaps
FUNC(mlkem_encaps)
    call    frame_enter                 // copies m
    push    rax                         // c, popped for c1
    push    rsi                         // K: G's output, consumed by hash_G
    push    rdx
    pop     rsi                         // ek (still in RDX)
    call    hash_H                      // m || H(ek), RDI = header+64
    call    hash_G                      // (K, r) = G(m || H(ek)): K goes to
                                        // the pushed output, RSI = r

    //  y[j] = CBD_eta1(PRF(r, j)), then e1[i] = CBD_eta2(PRF(r, K+i)) and
    //  e2 = CBD_eta2(PRF(r, 2K))
    mov     bh, 2*MLKEM_K+1             // y[K], e1[K], e2 from RBP
    call    cbd_vec                     // RSI = r; RDI ends after e2

    //  Rows from the last accumulator down, as in keygen: first
    //  v = e2 + t*y (BL = K: dot decodes t from ek and leaves R10 at rho),
    //  then u[i] = e1[i] + A^T*y for BL = i = K-1..0; dot runs BH = j and
    //  EBX is the seed suffix rho || i || j.
    mov     bl, MLKEM_K
.Lenc_row:
    call    dot                         // returns BH = K, RSI = accumulator
    dec     bl
    jns     .Lenc_row

    //  c1[i] = ByteEncode_du(Compress_du(u[i])), i = 0..K-1 (BH)
    pop     rdi                         // c
.Lenc_c1:
    mov     cl, DU                      // RCX = 0 after dot, then DU
    call    poly_compress_encode        // RSI ends at the next polynomial
    dec     bh
    jnz     .Lenc_c1
    //  v += Decompress_1(m), m at RBP-HDR. The row loop left BL = 0xFF
    //  and the c1 loop BH = 0 (poly_compress_encode keeps RBX), so EBX =
    //  255..0 indexes the coefficients.
.Lenc_mu:
    bt      dword ptr [rbp-HDR], ebx
    jnc     .Lenc_mu_next
    add     word ptr [rsi+rbx*2], (Q+1)/2
.Lenc_mu_next:
    dec     ebx
    jns     .Lenc_mu
    mov     cl, DV                      // ECX = DV (it was DU)
    //  Fall through to poly_compress_encode, which returns to frame_enter.
ENDF(mlkem_encaps)

//  poly_compress_encode: ByteEncode_d(Compress_d) of the polynomial at RSI
//  into RDI (32*d bytes), RCX = d (all 64 bits: LOOP counts the bits of
//  each coefficient from it). poly_encode is the same entry with RCX = 0
//  (keygen, after ntt): ByteEncode_12 without compression (d = 12 is
//  loaded per coefficient; the pushed RCX = 0 is restored). Compression
//  needs R8D = q (left by sample_cbd's exit in dot, in encaps and decaps;
//  nothing after dot changes R8). Inputs may exceed q (up to q + 1664 for
//  v + mu; adding q adds 2^d, which the packer drops). Compress_d(x) =
//  round(x * 2^d / q) is computed exactly as floor((x * 2^d + (q-1)/2) / q):
//  q is odd, so there are no ties. Only the low d bits of each value are
//  packed.
//  Each bit is rotated into the top of the output byte at RDI (RCR through
//  the carry: after 8 bits all of its old contents are gone); BL += 32
//  carries every 8 bits, then RDI moves to the next byte (SCASB). BH counts
//  256 coefficients. The output byte is read and rotated before it is
//  complete; in place (decaps' m', d = 1) byte k changes while
//  coefficients 8k..8k+7 are processed, and it holds part of coefficient
//  k/2 <= 8k, which has been loaded by then.
//  Returns RDI = end of output, RSI = end of the polynomial and RCX = d
//  (0 for poly_encode); keeps RBX, clobbers RAX and RDX.

FUNC(poly_compress_encode)
FUNC(poly_encode)
    push    rbx
    xor     ebx, ebx
.Lenc_coef:
    lodsw
    push    rcx
    jrcxz   .Lenc_raw                   // ByteEncode_12: bits 16+ unused
    cwde                                // x < 2^15
    shl     eax, cl                     // < 2^24
    add     eax, (Q-1)/2
    cdq
    div     r8d
.Lenc_bit:
    shr     eax, 1
    rcr     byte ptr [rdi], 1
    add     bl, 32
    jnc     .Lenc_next
    scasb                               // RDI += 1 (flags unused)
.Lenc_next:
    loop    .Lenc_bit
    pop     rcx
    dec     bh
    jnz     .Lenc_coef
.Lpe_ret:
    pop     rbx
    ret
.Lenc_raw:
    mov     cl, 12                      // RCX = 0 is pushed and restored
    jmp     .Lenc_bit
ENDF(poly_encode)
ENDF(poly_compress_encode)

//  mlkem_decaps: ML-KEM.Decaps_internal(dk, c) (Algorithm 18) with
//  K-PKE.Decrypt (Algorithm 15) inlined, in the normal domain. Re-encryption
//  calls mlkem_encaps, which recomputes h = H(ek) instead of reading it
//  from dk; the two are equal for every dk that passes the FIPS 203
//  decapsulation key check. Variable time.
//  RDI = shared key K, RSI = ciphertext, RDX = dk.

//  DEC_NEG * 65536 + 3303 = -3303 = -128^-1 (mod q): see invntt and
//  sample_cbd.
#define DEC_NEG     0x852

//  Frame: z in the header (frame_enter copies it), u[K] from RBP and the
//  accumulator, as in encaps; afterwards m' || c' from the accumulator
//  and, on a mismatch, z || c from the header.

    .globl  mlkem_decaps
FUNC(mlkem_decaps)
    lea     rcx, [rdx+768*MLKEM_K+64]   // z, after s^ || ek || H(ek)
    call    frame_enter
    push    rax                         // K: J's output, or popped at the end
    push    rdi                         // header+32: c goes after z
    push    rsi                         // c, for the compare and the copy
    push    rax                         // K, for encaps

    //  u[i] = Decompress_du(ByteDecode_du(c1[i])) at RBP, as y in encaps
    push    rbp
    pop     rdi                         // RSI = c (frame_enter keeps it)
    //  then v = Decompress_dv(ByteDecode_dv(c2)) into the accumulator, and
    //  poly_decode falls into dot. BH = K polynomials with DU before v;
    //  BL = 0 counts coefficients, then BL = K selects s from dk in dot.
    //  EBX bits 16+ = DEC_NEG make the invntt after sample_cbd's decode
    //  use c = -128^-1: dot adds (-s)^T u. (encaps, called below, saves
    //  and clears EBX.)
    mov     ebx, DEC_NEG << 16 | MLKEM_K << 8
    mov     cl, DU                      // RCX = 0 after frame_enter
    call    poly_decode                 // u[K], v, then dot: RSI = acc

    //  acc = v - s^T u = w; m' = ByteEncode_1(Compress_1(w)), encoded in
    //  place at the accumulator (see poly_compress_encode).
    mov     cl, 1                       // RCX = 0 after dot
    push    rdi                         // m'
    call    poly_compress_encode        // RDI = m' + 32, where c' goes

    //  (K, c') = Encaps_internal(ek, m') writes K straight to the output
    //  and c' after m'; on a mismatch K = J(z || c) (implicit rejection).
    pop     rcx                         // m', copied before c' is written
    pop     rsi                         // K
    mov     rdx, r10                    // ek inside dk
    push    rdi                         // c'
    call    mlkem_encaps                // returns RCX = DV (poly_encode)
    pop     rdi                         // c'
    pop     rsi
    push    rsi                         // c
    mov     cl, CT_BYTES/8
    repe    cmpsq                       // ZF = 1: equal; kept until the JE
    pop     rsi                         // c
    pop     rdi                         // header+32
    mov     cl, CT_BYTES/8              // REPE left RCX < 256
    rep     movsq                       // z || c in the frame (unused if equal)
    je      .Lpe_ret                    // equal: K = K' (encaps wrote it);
                                        // pops the K pointer into RBX
    mov     r8w, SHAKE256               // K = J(z || c): the copy ended
                                        // at RDI, the output K is at [RSP]
    //  Fall through to hash_hdr (Keccak section follows).
ENDF(mlkem_decaps)

//  ========================================================================
//  Keccak
//  ========================================================================

//  hash_hdr: hash the bytes from RBP-HDR up to RDI with mode R8W into the
//  output at [RSP], via xhash_out with the Keccak state at RDI, right
//  after the input (no return address of its own: mlkem_decaps falls in
//  for J, the state then lies in decaps' dead frame; hash_G jumps in).
//  hash_G: G = SHA3-512 of the bytes from RBP-HDR up to RDI, with the
//  state at RDI (its 200 bytes reach into the vector, which is sampled
//  later: s or y), so G's 64 bytes are in place there; the first 32 are
//  also copied to the output pointer that the caller pushed before the
//  CALL (consumed; keygen pushes RDI itself, encaps its K output).
//  Returns RSI = RDI + 32 (the second half: sigma or r).
//  hash_H: H = SHA3-256 (32 bytes) of the EK_BYTES-byte ek at RSI into
//  RDI, with the state at RBP (s, dead in keygen's tail, or y, sampled
//  later in encaps).
//  Both return RDI = output + 32, like xhash_out.

FUNC(hash_hdr)
    lea     rsi, [rbp-HDR]
    mov     edx, edi
    sub     edx, esi                    // input length
    //  Fall through to xhash_out.
ENDF(hash_hdr)

//  xhash_out: one-shot SHA3 / SHAKE; the first 32 bytes of the first
//  block are copied to the output pointer at [RSP] (popped on return).
//  RDI = Keccak state (200 bytes of scratch, not overlapping the input;
//  the output may be the state itself), RSI = input, RDX = input length,
//  R8W = mode, RCX < 256. Returns RSI = state + 32, RDI = output + 32
//  and RCX = 0; preserves R8-R10.

FUNC(xhash_out)
    call    keccak_absorb               // RCX = 0
    pop     rsi                         // output
    xchg    rsi, rdi                    // RSI = state, RDI = output
    mov     cl, 32                      // RCX = 0 after Keccak
    rep     movsb
    ret
ENDF(xhash_out)

FUNC(hash_G)
    mov     r8w, SHA3_512
    pop     rax
    pop     rdx                         // the caller's output pointer
    push    rax
    push    rdx                         // on top, the return address below
    jmp     hash_hdr
ENDF(hash_G)

FUNC(hash_H)
    mov     r8w, SHA3_256
    mov     edx, EK_BYTES
    push    rdi                         // output
    push    rbp
    pop     rdi                         // state
    jmp     xhash_out
ENDF(hash_H)


//  keccak_absorb: RDI = state (200 bytes, gets zeroed), RSI = input,
//  RDX = input length, R8W = mode (upper bits ignored), RCX < 256.
//  Absorbs, pads and permutes, so the first output block is ready in the
//  state. Preserves RDI, R8, R9, R10.

FUNC(keccak_absorb)
    push    rdi
    mov     cl, 200                     // RCX < 256 in all callers
    xor     eax, eax
    rep     stosb                       // 200 bytes
    pop     rdi                         // RCX = 0: position in block
.Lab_byte:
    dec     edx
    js      .Lab_pad
    lodsb
    xor     [rdi+rcx], al
    inc     ecx
    cmp     cl, r8b                     // block full?
    jne     .Lab_byte
    push    rsi
    push    rdx
    call    keccak_f1600
    pop     rdx
    pop     rsi                         // RCX = 0 after Keccak
    jmp     .Lab_byte
.Lab_pad:
    mov     eax, r8d
    xor     [rdi+rcx], ah               // domain + first padding bit
    mov     cl, al                      // RCX = rate (RCX < 256 already)
    xor     byte ptr [rdi+rcx-1], 0x80  // last padding bit
    //  Fall through to the permutation.
ENDF(keccak_absorb)

//  keccak_f1600: permute the 25-lane state at RDI (lane x + 5*y).
//  Round constants come from the LFSR x^8+x^6+x^5+x^4+1 in DH. Its
//  reciprocal right-shift form starts at 0x1d and uses 0x8e feedback: the
//  carry from SHR is the next iota bit. Each round builds its scratch
//  values with pushes; PUSH RBP / LEAVE frame them (about 500 bytes of
//  stack). Needs RCX < 256 (all callers have a small count there).
//  Clobbers RAX, DL, DH, RSI, R11 and flags (EDX bits 16+ are kept);
//  returns RCX = 0 (xhash_out, keccak_absorb and sample_cbd rely on it),
//  DL = 0 (the CBD coefficient count), RSI = RDI + 200 and EAX = 127
//  (SampleNTT relies on both: get_bits keeps EAX bits 12+, which must be
//  0; EAX bits 16+ = 0 in CBD). Preserves RDI, R8-R10 and all
//  callee-saved registers.

FUNC(keccak_f1600)
    mov     dh, 0x1d
.Lkf_round:
    push    rbp
    push    rsp
    pop     rbp
    //  theta: C[x] by a 5-deep running xor over lanes 24..0 (lanes
    //  24..20 start the columns, so they skip the xor)
    mov     cl, 25                      // RCX < 256 on entry, 0 later
.Lkf_col:
    mov     rax, [rdi+rcx*8-8]
    cmp     cl, 21
    jae     .Lkf_col_top
    xor     rax, [rsp+32]
.Lkf_col_top:
    push    rax
    loop    .Lkf_col
    //  [rsp + 8x] = C[x]. Lanes 24..0: D = C[x-1] ^ rol(C[x+1], 1),
    //  re-pushing C[x] keeps the window C[x+1] at [rsp], C[x-1] at +24.
    //  RAX = [rsp] on every entry (the last push, C[0] at first).
    mov     cl, 25
.Lkf_theta:
    rol     rax, 1
    xor     rax, [rsp+24]
    xor     [rdi+rcx*8-8], rax
    mov     rax, [rsp+32]
    push    rax
    loop    .Lkf_theta
    //  rho and pi: lane i = x + 5y moves to y + 5((2x + 3y) mod 5). EAX =
    //  2i and ESI = 2y, so EAX - ESI = 2x + 8y = 2x + 3y mod 5. Start one
    //  step early at (x, y) = (1, 1), whose successor is lane 1: that step
    //  stores junk there with rotation 0 and picks up A[1]; the last step
    //  overwrites the junk. CL = (t+1)t/2 rotation, CH = step count.
    mov     dl, 5                       // divisor here, chi's row count later
    push    2
    pop     rax                         // 2i = 2 for lane 6 (= 2x mod 5)
    push    rax
    pop     rsi                         // 2y = 2
.Lkf_rho:
    rol     r11, cl
    sub     eax, esi                    // 2x + 8y < 256, = 2x + 3y mod 5
    div     dl                          // AH = y' = (2x + 3y) mod 5, AL < 128
    shr     eax, 7                      // 2y'
    xchg    eax, esi                    // ESI = 2y', EAX = 2y = 2x'
    lea     eax, [rax+rsi*4]
    add     eax, esi                    // 2i' = 2x' + 10y'
    xchg    [rdi+rax*4], r11
    inc     ch
    add     cl, ch
    cmp     ch, 25
    jne     .Lkf_rho
    //  chi: per row push A1 A0 A4 A3 A2 A1, then pop A[x+1] and read A[x+2]
    push    rdi
    pop     rsi                         // DL = 5 rows (set before rho)
.Lkf_chi_row:
    push    qword ptr [rsi+8]
    push    qword ptr [rsi]
    push    4
    pop     rcx
.Lkf_chi_push:
    push    qword ptr [rsi+rcx*8]
    loop    .Lkf_chi_push
    mov     cl, 5
.Lkf_chi:
    pop     rax
    not     rax
    and     rax, [rsp]
    xor     [rsi], rax
    lodsq
    loop    .Lkf_chi
    dec     dl
    jnz     .Lkf_chi_row
    //  iota: LFSR output bit j toggles lane bit 2^j - 1, j = 0..6
    xor     eax, eax                    // bit position 0, 1, 3, .., 63
    mov     cl, 7
.Lkf_iota:
    shr     dh, 1                       // emit old bit 0 into carry
    jnc     .Lkf_iota_skip
    xor     dh, 0x8e                    // reciprocal feedback if bit was 1
    btc     dword ptr [rdi], eax        // dword form reaches bit 63 too
.Lkf_iota_skip:
    stc
    adc     eax, eax                    // 2 * EAX + 1
    loop    .Lkf_iota
    leave
    cmp     dh, 0xe2                    // LFSR state after 24 rounds
    jne     .Lkf_round
    ret
ENDF(keccak_f1600)

    .section .note.GNU-stack, "", @progbits
