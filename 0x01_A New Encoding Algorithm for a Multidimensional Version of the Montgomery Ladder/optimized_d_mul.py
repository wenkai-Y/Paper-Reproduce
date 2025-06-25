# the simultaneous differential multidimensional scalar point multiplication algorithm
# d-MUL 同步微分多维标量点乘算法
# from numpy.random.c_distributions import random_exponential
import math
from traceback import print_tb

# from numpy.random.c_distributions import random_standard_normal
from scapy.layers.msrpce.ept import prot_and_addr_t
from scapy.layers.spnego import NEGOEX_BYTE_VECTOR
import time
import timeit
from No_02_scalar_multiplication import point_add_jacobian
from No_02_scalar_multiplication import point_double_jacobian
from No_02_scalar_multiplication import montgomery_ladder_jacobian
import random
G_X = "32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7"
G_x = int(G_X, 16)
G_Y = "BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0"
G_y = int(G_Y, 16)
p = "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF"
p = int(p, 16)
n = "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123"
n = int(n, 16)

def Sanitize(d, a, P):
    for i in range(d):
        if a[i] < 0:
            a[i] = -a[i]
            P[i][1] = -P[i][1]
    return a, P

def ChooseSeq(d, a):
    evens = []
    odds  = []
    for i in range(d):
        if a[i] % 2 == 0:
            evens.append(i+1)
        else:
            odds.append(i+1)
    # 这里就不乱序了哈
    return odds + evens

def from_jacobian_to_affine(Q):
    if Q != (1, 1, 0):
        X, Y, Z = Q
        # 将结果从雅可比坐标转换回标准坐标
        X_res = X * (pow(pow(Z, 2, p), -1, p)) % p
        Y_res = Y * (pow(pow(Z, 3, p), -1, p)) % p
        return (X_res, Y_res)
    else:
        print("结果是无穷远点")
        return None


# P是雅可比坐标
def Optimized_d_MUL(d, _a, _P):
    # step-1
    a, P = Sanitize(d, _a, _P)
    # a = _a
    # step-2

    sigma = ChooseSeq(d, a)
    # print("sigma: ", sigma)

    # step-3
    ell = max((a[i].bit_length() for i in range(d)))
    ahat = [a[i] + (a[i]%2) - 1 for i in range(d)]
    # print("ahat: ", ahat)
    b = []
    for i in range(d):
        bits = list(map(int, bin(ahat[i])[2:]))
        bits = [0] + [0] * (ell - len(bits)) + bits
        b.append(bits)
    # print("b: ", b)

    # step-4
    r = [0] * (ell * d)
    for k in range(ell, 0, -1):
        L0, L1 = [], []
        for i in range(1, d+1, 1):
            idx = (k - 1) * d + i - 1
            b1 = b[sigma[i - 1]-1][k - 1]
            b2 = b[sigma[i - 1]-1][k]
            r[idx] = b1 ^ b2
            if b1 ^ b2 == 0:
                L0.append(sigma[i - 1])
            else:
                L1.append(sigma[i - 1])
        sigma = list(reversed(L1)) + L0
    # print("r: ", r)
    # print("sigma: ", sigma)

    # step-5
    Q = [(1, 1, 0) for _ in range(d + 1)]
    R = [(1, 1, 0) for _ in range(d + 1)]
    for i in range(1, d+1, 1):
        # 点加
        Q[i] = point_add_jacobian(Q[i-1], P[sigma[i-1]-1])

    # step-6
    for k in range(1, ell+1, 1):
        # 找到目标字符串的汉明重量
        h = 1
        for index in range((k-1)*d+1, (k-1)*d+1+d, 1):
            h += r[index-1]
        x, y = h, h
        # 一次倍加
        R[0] = point_double_jacobian(Q[h-1])
        for i in range(1, d+1, 1):
            # d次点加
            x -= r[(k-1)*d+i-1]
            y += 1 - r[(k-1)*d+i-1]
            R[i] = point_add_jacobian(Q[x-1], Q[y-1])
        # 切记要深拷贝 :(
        Q = R.copy()

    # step-7
    h = 1
    for i in range(d):
        h += a[i] % 2
    # print("Q: ", Q)
    # print("h: ", h)
    return Q[h-1]

def easy_d_MUL(d, _a, _P):
    Q = []
    R = (1, 1, 0)
    for i in range(d):
        Q.append(montgomery_ladder_jacobian(_a[i], _P[i]))
    for i in range(d):
        R = point_add_jacobian(R, Q[i])
    return R











if __name__ == '__main__':

    G = (G_x, G_y, 1)
    a = int('0x9EBBBE0681EB4C77978CD7FD2A59935991C2BB595437B52FD5A108756F706E61', 16)
    b = int('0xFC57E5325D877E5E72C3E59FE8DA572075A22D259C94FD0460BD0B2AA1483498', 16)
    P1 = (a, b, 1)

    for j in range(1000):
        print("No.", j)

        a = [random.randint(1, n - 1) for _ in range(4)]
        # P= [G, G, G, G]
        P= [G, P1, G, P1]
        Q = Optimized_d_MUL(4, a, P)
        # Q = Optimized_d_MUL(4, [1, 4, 5, 6], [10, 2, 45, 9])
        Q_easy = easy_d_MUL(4, a, P)
        Q = from_jacobian_to_affine(Q)
        Q_easy = from_jacobian_to_affine(Q_easy)
        if (Q != Q_easy):
            print("error")
        else:
            print("success")

    # a = [random.randint(1, n - 1) for _ in range(10)]
    # P = [G, P1, G, P1, G, P1, G, P1, G, P1]
    #
    # time_taken = timeit.timeit("Optimized_d_MUL(10, a, P)", globals=globals(), number=1000)
    # print("平均执行时间：", time_taken, "毫秒")
    # time_taken = timeit.timeit("easy_d_MUL(10, a, P)", globals=globals(), number=1000)
    # print("平均执行时间：", time_taken, "毫秒")
