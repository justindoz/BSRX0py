import struct
from copy import deepcopy

def U8TO32(p):
    """Converts 4 bytes to a 32-bit unsigned integer (big-endian)."""
    return (p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3]

def U32TO8(p, v):
    """Converts a 32-bit unsigned integer to 4 bytes (big-endian)."""
    p[0] = (v >> 24) & 0xFF
    p[1] = (v >> 16) & 0xFF
    p[2] = (v >> 8) & 0xFF
    p[3] = v & 0xFF

def ROTR32(x, n):
    """Performs a 32-bit right rotation."""
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

# Constants as defined in the original code
d_blake_sigma = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
    [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
    [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
    [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
    [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
    [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
    [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8]
]

d_blake_cst = [
    0x243F6A88, 0x85A308D3, 0x13198A2E, 0x03707344,
    0xA4093822, 0x299F31D0, 0x082EFA98, 0xEC4E6C89,
    0x452821E6, 0x38D01377, 0xBE5466CF, 0x34E90C6C,
    0xC0AC29B7, 0xC97C50DD, 0x3F84D5B5, 0xB5470917
]

class BlakeState:
    """Represents the state of the BLAKE hash function."""
    def __init__(self):
        self.h = [0] * 8      # Hash state
        self.s = [0] * 4      # Salt
        self.t = [0] * 2      # Counters
        self.buflen = 0       # Buffer length in bits
        self.nullt = 0        # Null flag
        self.buf = bytearray(64)  # Data buffer

def BLAKE_G(v, m, sigma, cst, i, a, b, c, d, e):
    """Implements the G mixing function."""
    v[a] = (v[a] + ((m[sigma[i][e]] ^ cst[sigma[i][e+1]])) + v[b]) & 0xFFFFFFFF
    v[d] = ROTR32(v[d] ^ v[a], 16)
    v[c] = (v[c] + v[d]) & 0xFFFFFFFF
    v[b] = ROTR32(v[b] ^ v[c], 12)
    v[a] = (v[a] + ((m[sigma[i][e+1]] ^ cst[sigma[i][e]])) + v[b]) & 0xFFFFFFFF
    v[d] = ROTR32(v[d] ^ v[a], 8)
    v[c] = (v[c] + v[d]) & 0xFFFFFFFF
    v[b] = ROTR32(v[b] ^ v[c], 7)

def cn_blake_compress(S, block):
    """Compression function for BLAKE."""
    v = [0] * 16
    m = [0] * 16

    # Convert block bytes to 32-bit words
    for i in range(16):
        m[i] = U8TO32(block[i*4:(i+1)*4])

    # Initialize the state vector
    for i in range(8):
        v[i] = S.h[i]
    v[8]  = S.s[0] ^ 0x243F6A88
    v[9]  = S.s[1] ^ 0x85A308D3
    v[10] = S.s[2] ^ 0x13198A2E
    v[11] = S.s[3] ^ 0x03707344
    v[12] = 0xA4093822
    v[13] = 0x299F31D0
    v[14] = 0x082EFA98
    v[15] = 0xEC4E6C89

    if S.nullt == 0:
        v[12] ^= S.t[0]
        v[13] ^= S.t[0]
        v[14] ^= S.t[1]
        v[15] ^= S.t[1]

    # 14 Rounds
    for i in range(14):
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 0, 4, 8, 12, 0)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 1, 5, 9, 13, 2)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 2, 6, 10, 14, 4)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 3, 7, 11, 15, 6)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 3, 4, 9, 14, 14)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 2, 7, 8, 13, 12)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 0, 5, 10, 15, 8)
        BLAKE_G(v, m, d_blake_sigma, d_blake_cst, i, 1, 6, 11, 12, 10)

    # Update the hash state
    for i in range(16):
        S.h[i % 8] ^= v[i]
    for i in range(8):
        S.h[i] ^= S.s[i % 4]

def cn_blake_update(S, data, datalen_bits):
    """Updates the BLAKE state with new data."""
    left = S.buflen >> 3  # Number of bytes in buffer
    fill = 64 - left      # Remaining space in buffer

    data_len_bytes = datalen_bits >> 3

    if left and ((data_len_bytes & 0x3F) >= fill):
        # Fill the buffer and compress
        S.buf[left:left+fill] = data[:fill]
        S.t[0] = (S.t[0] + 512) & 0xFFFFFFFF
        if S.t[0] == 0:
            S.t[1] = (S.t[1] + 1) & 0xFFFFFFFF
        cn_blake_compress(S, S.buf)
        data = data[fill:]
        datalen_bits -= (fill << 3)
        left = 0

    # Compress blocks of 64 bytes (512 bits)
    while datalen_bits >= 512:
        S.t[0] = (S.t[0] + 512) & 0xFFFFFFFF
        if S.t[0] == 0:
            S.t[1] = (S.t[1] + 1) & 0xFFFFFFFF
        cn_blake_compress(S, data[:64])
        data = data[64:]
        datalen_bits -= 512

    # Buffer remaining data
    if datalen_bits > 0:
        num_bytes = datalen_bits >> 3
        S.buf[left:left+num_bytes] = data[:num_bytes]
        S.buflen = (left << 3) + datalen_bits
    else:
        S.buflen = 0

def cn_blake_final(S, digest):
    """Finalizes the BLAKE hash computation."""
    padding = bytes([
        0x80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    ])

    pa = bytes([0x81])
    pb = bytes([0x01])
    msglen = bytearray(8)

    lo = (S.t[0] + S.buflen) & 0xFFFFFFFF
    hi = S.t[1]
    if lo < S.buflen:
        hi = (hi + 1) & 0xFFFFFFFF

    U32TO8(msglen, hi)
    U32TO8(msglen, lo)  # Overwrites first 4 bytes
    U32TO8(msglen[4:], hi)
    U32TO8(msglen[4:], lo)

    if S.buflen == 440:
        S.t[0] = (S.t[0] - 8) & 0xFFFFFFFF
        cn_blake_update(S, pa, 8)
    else:
        if S.buflen < 440:
            if S.buflen == 0:
                S.nullt = 1
            S.t[0] = (S.t[0] - (440 - S.buflen)) & 0xFFFFFFFF
            cn_blake_update(S, padding[:(440 - S.buflen)], (440 - S.buflen) * 8)
        else:
            S.t[0] = (S.t[0] - (512 - S.buflen)) & 0xFFFFFFFF
            cn_blake_update(S, padding[:(512 - S.buflen)], (512 - S.buflen) * 8)
            S.t[0] = (S.t[0] - 440) & 0xFFFFFFFF
            cn_blake_update(S, padding[1:441], 440 * 8)
            S.nullt = 1
        cn_blake_update(S, pb, 8)
        S.t[0] = (S.t[0] - 8) & 0xFFFFFFFF

    S.t[0] = (S.t[0] - 64) & 0xFFFFFFFF
    cn_blake_update(S, msglen, 64)

    # Convert the final hash state to bytes
    for i in range(8):
        U32TO8(digest[i*4:(i+1)*4], S.h[i])

def cn_blake(in_bytes, inlen_bits):
    """Main BLAKE hash function."""
    S = BlakeState()

    # Initialize hash state
    S.h[0] = 0x6A09E667
    S.h[1] = 0xBB67AE85
    S.h[2] = 0x3C6EF372
    S.h[3] = 0xA54FF53A
    S.h[4] = 0x510E527F
    S.h[5] = 0x9B05688C
    S.h[6] = 0x1F83D9AB
    S.h[7] = 0x5BE0CD19

    # Initialize other state variables
    S.t = [0, 0]
    S.buflen = 0
    S.nullt = 0
    S.s = [0, 0, 0, 0]

    # Update state with input data
    cn_blake_update(S, in_bytes, inlen_bits)

    # Finalize and get the digest
    digest = bytearray(32)
    cn_blake_final(S, digest)

    return bytes(digest)

# Example usage
if __name__ == "__main__":
    # Example input
    input_data = b"Hello, BLAKE!"
    input_length_bits = len(input_data) * 8

    # Compute BLAKE hash
    hash_digest = cn_blake(input_data, input_length_bits)

    # Print the hash in hexadecimal format
    print("BLAKE Hash:", hash_digest.hex())
