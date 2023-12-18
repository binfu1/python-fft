#! /usr/bin/python3
# -*- encoding:utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, cos
import os
import math
import cmath

#对原始信号进行采样
fs = 32  #采样频率， 要大于信号频率的两倍
t = np.arange(0, 1, 1.0/fs)  #1秒采样fs个点
N = len(t)
freq = np.arange(N) #频率counter
f = 3*cos(2*pi*2*t) #离散化后的f[n]，角频率为2，周期为1/2=0.5

# 直接调用函数库进行计算
fft = np.fft.fft(f) #快速傅里叶变换
# 换算成实际的振幅
#fft = fft/(N/2)  
#fft[0] = fft[0]/2
for i in range(0, N):
  if (abs(fft[i]) > 0.1):
    print('频率%d处的幅度为%.1f' %(i, abs(fft[i])))

# 根据DFT公式原理实现的DFT计算
# W = np.zeros((N,N), complex)
dft = np.zeros(N, complex)
for k in range(0, N):
  for n in range(0, N):
    #dft[k] = dft[k] + (1/N)*f[n]*np.exp(-2j*pi*k*n/N)
    #W[k,n] = np.exp(-2j*pi*k*n/N)
    #dft[k] = dft[k] + f[n]*W[k,n]
    dft[k] = dft[k] + f[n]*np.exp(-2j*pi*k*n/N)
  #if (abs(dft[k]) > 0.1):
  #  print('频率%d处的幅度为%.1f' %(k, abs(dft[k])))
'''
filename = 'W'
os.remove(filename)
with open(filename, 'w') as file:
  for k in range(N):
    for n in range(N):
      if n==N-1:
        file.write(str(W[k,n])+'\n\n')
      elif n%4==3:
        file.write(str(W[k,n])+'\n')
      else:
        file.write(str(W[k,n])+' ')
file.close
'''

# 根据FFT公式原理实现的FFT计算
bit =int(math.log2(N))
# 通过位反转重新排列输入数据
f1 = np.zeros(N, complex)
for i in range(N):
  ir = int(f"{i:0{bit}b}"[::-1], 2)
  f1[i] = f[ir]
for s in range(1, N.bit_length()):
  m = 2 ** s
  for k in range(0, N, m):
    for j in range(0, m // 2):
      if (k + j + m // 2==N):
        break
      W= cmath.exp(-2j * cmath.pi * j*(2 ** (N.bit_length()-s-1)) / N)
      x2 = W* f1[k + j + m // 2]
      x1 = f1[k + j]
      f1[k + j] = x1 + x2
      f1[k + j + m // 2] = x1 - x2

# 根据IFFT公式原理实现的IFFT计算
f2 = f1
for i in range(N): #频域共轭
  f2_real =  f2[i].real
  f2_imag = -f2[i].imag
  f2[i] = complex(f2_real, f2_imag)
bit =int(math.log2(N))
# 通过位反转重新排列输入数据
f3 = np.zeros(N, complex)
f4 = np.zeros(N)
for i in range(N):
  ir = int(f"{i:0{bit}b}"[::-1], 2)
  f3[i] = f2[ir]
for s in range(1, N.bit_length()):
  m = 2 ** s
  for k in range(0, N, m):
    for j in range(0, m // 2):
      if (k + j + m // 2==N):
        break
      W= cmath.exp(-2j * cmath.pi * j*(2 ** (N.bit_length()-s-1)) / N)
      x2 = W* f3[k + j + m // 2]
      x1 = f3[k + j]
      f3[k + j] = x1 + x2
      f3[k + j + m // 2] = x1 - x2
for i in range(N): #时域共轭
  f3_real = f3[i].real/N
  f3_imag = -f3[i].imag/N
  f3[i] = complex(f3_real, f3_imag)
  f4[i] = f3[i].real

# 绘制时域图像
plt.figure(figsize=(10,6))
plt.subplots_adjust(wspace=0.3, hspace=0.5)
plt.subplot(2,2,1)
plt.plot(t, f, 'bo-')
plt.title('signal')
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.subplot(2,2,2)
plt.plot(freq, abs(fft), 'ro-')
plt.title('fft-lib')
plt.xlabel('freq (Hz)')
plt.ylabel('amplitude')
plt.subplot(2,2,3)
#plt.plot(freq, abs(dft), 'go-')
#plt.title('dft')
#plt.xlabel('freq (Hz)')
#plt.ylabel('amplitude')
plt.plot(t, f4, 'go-')
plt.title('ifft')
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.subplot(2,2,4)
plt.plot(freq, abs(f1), 'mo-')
plt.title('fft-alg')
plt.xlabel('freq (Hz)')
plt.ylabel('amplitude')
plt.show()