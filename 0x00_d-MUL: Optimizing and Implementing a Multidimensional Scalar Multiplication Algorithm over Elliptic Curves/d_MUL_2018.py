# 对2018年的这篇d-MUL论文的算法进行复现
# d-MUL: Optimizing and Implementing a Multidimensional Scalar Multiplication Algorithm over Elliptic Curves

import random
from typing import List
from No_02_scalar_multiplication import point_add_jacobian
from No_02_scalar_multiplication import point_double_jacobian
from No_02_scalar_multiplication import montgomery_ladder_jacobian

p = "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF"
p = int(p, 16)
IDE = (1, 1, 0)

# Algorithm 1: d-MUL scalars
def d_MUL_scalars(l: int, r: List[int], tau: List[int], v: List[int]) -> List[int]:
    """
    实现 d-MUL scalar 生成算法。

    参数:
        l   : 标量比特长度（每个标量的位数）
        r   : 比特串（长度为 l * d）
        tau : 一个对 {0,...,d-1} 的排列（长度为 d）
        v   : 长度为 d 的比特串，用于修正最终结果

    返回:
        a  : 长度为 d 的整数列表 [a1, ..., ad]
    """
    d = len(tau)
    B = [[0] * d]
    for i in range(d):
        new_row = B[i][:]
        new_row[tau[i]] += 1
        B.append(new_row)
    for i in range(l):
        segment = r[i*d : (i+1)*d]
        h = sum(segment)
        x, y = h, h
        A = []
        A.append([2 * val for val in B[h]])
        for j in range(d):
            if segment[j] == 1:
                x -= 1
            else:
                y += 1
            A.append([B[x][k] + B[y][k] for k in range(d)])
        B = A

    return [B[d][i] - v[i] for i in range(d)]


# Algorithm 2: simplified_d_MUL
def simplified_d_MUL(l, P, r, tau, v):
    """
    实现简化版 d-MUL 多维标量乘法算法（Algorithm 2）

    参数:
        l   : 标量的比特长度（标量每位的长度）
        P   : 椭圆曲线点列表 [P1, ..., Pd]每个 Pi 是椭圆曲线上的点【雅可比坐标】
        r   : 比特串（长度为 l * d），用于控制加法链结构
        tau : 一个对 {0, ..., d-1} 的排列，用于构造初始状态点序列
        v   : 长度为 d 的比特串，表示输出结果中应被减去的“偏移点”

    返回:
        a   : 对应 Q 的 scalars 向量 [a1, ..., ad]，由 Algorithm 1 得出
        T   : 椭圆曲线上的点 T = a1*P1 + ... + ad*Pd 实际上是由结构性加法链构造得到的 Q[d]，再减去 v 对应的点
    """
    d = len(P)
    Q = [IDE]
    for i in range(d):
        Q.append(point_add_jacobian(Q[i], P[tau[i]]))
    for i in range(l):
        segment = r[i * d : (i + 1) * d]
        h = sum(segment)
        x, y, = h, h
        R = [point_double_jacobian(Q[h])]
        for j in range(d):
            if segment[j] == 1:
                x -= 1
            else:
                y += 1
            R.append(point_add_jacobian(Q[x], Q[y]))
        Q = R
    T = Q[d]
    for i in range(d):
        if v[i] == 1:
            T = point_add_jacobian(T, (P[i][0], -P[i][1], P[i][2]))
    a = d_MUL_scalars(l, r, tau, v)
    return a, T

# BoolToInt
def BoolToInt(b):
    return 1 if b else 0
# Select
def Select(a, b, cond):
    return b if cond else a
# Xnor 同或门
def Xnor(a, b):
    return int(not(a ^ b))

#  Algorithm 4: d-MUL scalars (Optimized)
def d_MUL_scalars_Optimized(r, tau, l, d):
    k = [0 for _ in range(d)]
    for i in range(d):
        k[i] = 0
        delta = 1
        index = i
        for j in range(0, d * (l - 2), d):
            h = sum(r[j+t] for t in range(d))
            z = BoolToInt(h > index)
            k[i] = 2 * (k[i] + delta)
            delta = (1 - 2 * z) * delta
            q = index + 1 - h
            a = 0
            index = -1
            q = Select(q, -q, BoolToInt(q > 0)) + z
            for t in range(d):
                a = a + Xnor(z, r[j + t])
                index = Select(t, index, BoolToInt((a == q) and (index == -1)))
        idx = (l - 1) * d + tau[i] - 1
        idx = max(0, min(idx, len(r)-1))
        k[tau[i]] = 2 * k[i] + delta - r[idx]
    return k




# main
if __name__ == '__main__':
    l = 256  # 每个标量的 bit 长度
    d = 2  # 维度（标量个数）

    r = [random.randint(0, 1) for _ in range(l * d)]  # 随机比特串
    tau = random.sample(range(d), d)  # 随机排列（不重复）
    v = [random.randint(0, 1) for _ in range(d)]  # 随机偏移向量

    # print("输入参数：")
    # print("比特串 r:")
    # print(r)
    # print("排列 tau:", tau)
    # print("偏移 v   :", v)

    K = d_MUL_scalars_Optimized(r, tau, l, d)
    print("d_MUL_scalars_Optimized : 输出标量向量 [a1, a2, ..., ad]:")
    print(K)

    P1 = (int('0x8EE77087009544506DF55BDFE24112A5D4E8F556A0AC6876565EA45983DA3D72', 16),
          int('0xFB9B3F0C20270464E2B21531C8ADCD9AD981CDBD98909A3D1D1A9B3360CD1B11', 16),
          1)
    P2 = (int('0x3739C83AF9BC00D843803140011FB035D62242C635BA4F7DD394C636BBB6A874', 16),
          int('0xFBF5B20C2B6A7746BB9EC8DAA0F01E27B88641525429C3E35E0FAD54AEE76031', 16),
          1)
    P = [P1, P2]

    scalars, T = simplified_d_MUL(l, P, r, tau, v)

    print("输出标量向量 [a1, a2, ..., ad]:")
    print(scalars)
    if T != (1, 1, 0):
        X, Y, Z = T
        # 将结果从雅可比坐标转换回标准坐标
        X_res = X * (pow(pow(Z, 2, p), -1, p)) % p
        Y_res = Y * (pow(pow(Z, 3, p), -1, p)) % p
        print(f"kP_x = {hex(X_res)}")
        print(f"kP_y = {hex(Y_res)}")
    else:
        print("结果是无穷远点")

    # 测试即可
    Q1 = montgomery_ladder_jacobian(scalars[0], P[0])
    Q2 = montgomery_ladder_jacobian(scalars[1], P[1])
    Q = point_add_jacobian(Q1, Q2)
    if Q != (1, 1, 0):
        X, Y, Z = Q
        # 将结果从雅可比坐标转换回标准坐标
        X_res = X * (pow(pow(Z, 2, p), -1, p)) % p
        Y_res = Y * (pow(pow(Z, 3, p), -1, p)) % p
        print(f"kP_x = {hex(X_res)}")
        print(f"kP_y = {hex(Y_res)}")
    else:
        print("结果是无穷远点")
