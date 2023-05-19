import matplotlib.pyplot as plt

data_1 =  [6.594000, 24.799000, 96.153000, 382.254000, 852.564000, 1500.099000 ,2353.092000]
size =  [1000,     2000,      4000,      8000,       12000,      16000,       20000]

data_2 = [140.979999, 538.265745, 2216.818424,10335.557138, 21970.693473, 46447.083816, 67933.167603]
data_4 = [166.419392, 593.488430, 2343.432016, 10381.217610, 22821.818073, 43901.514795, 67110.991779]

plt.plot(size, data_1, label = "Non parallel")
plt.plot(size, data_2, label = "2 - threads")
plt.plot(size, data_4, label = "4 - threads")

plt.legend()
plt.grid()

plt.show()