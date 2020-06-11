import math

# iterative cooley tukey fft for dyadics
def fft(v):
    n, h  = len(v), len(v) >> 1
    old, new = v[:], [0] * n
    sublen, stride = 1, n

    while sublen < n:
        stride >>= 1
        for i in range(stride):
            for k in range(0, n, 2*stride):
                omega = math.exp(math.pi * k / n)
                new[i+(k>>1)]   = old[i+k] + omega * old[i+k+stride]
                
                new[i+(k>>1)+h] = old[i+k] - omega * old[i+k+stride]
        old, new = new, old
        sublen <<= 1

    return old
