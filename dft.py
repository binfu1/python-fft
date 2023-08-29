#! /usr/bin/python3
# -*- encoding:utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.pylab import mpl
import numpy as np
from numpy import pi, cos, sin

# 对原始信号进行采样
fs = 256  # 采样频率， 要大于信号频率的两倍
t = np.arange(0, 1, 1.0 / fs)  # 1秒采样fs个点
N = len(t)
freq = np.arange(N)  # 频率counter
f = 2 + 3 * cos(2 * pi * 50 * t) + 1.5 * cos(2 * pi * 75 * t)  # 离散化后的f[n]

# 直接调用函数库进行计算
F = np.fft.fft(f)  # 离散傅里叶变换
for i in range(0, N):
    if (abs(F[i]) > 0.25):
        print('%d is %d', i, abs(F[i]))

# 根据DFT公式原理，实现的DFT计算，做了/N的标准化
F2 = np.zeros(N, dtype=complex)  # F[n]
for k in range(0, N):  # 0,1,2,...,N-1
    for n in range(0, N):  # 0,1,2,...,N-1
        F2[k] = F2[k] + (1 / N) * f[n] * np.exp(-2j * pi * k * n / N)

# 绘制时域图像
fig, ax = plt.subplots(3, 1, figsize=(12, 12))
ax[0].plot(t, f, label='original signal')
ax[0].set_xlabel('time (s)')
ax[0].set_ylabel('amplitude')
ax[1].plot(freq, abs(F), 'r', label='use python library')
ax[1].set_xlabel('freq (Hz)')
ax[1].set_ylabel('amplitude')
ax[1].legend()
ax[2].plot(freq, abs(F2), 'r', label='compute by formula')
ax[2].set_xlabel('freq (Hz)')
ax[2].set_ylabel('amplitude')
ax[2].legend()
plt.show()
