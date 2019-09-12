import matplotlib.pyplot as plt
import math

residue, time, var1, var2 = [], [], [], []

with open("residueVsTime") as f:
    line = f.readlines()
    for l in line:
        l = l.split(",")
        l = [float(e) for e in l]
        residue.append(l[0])
        time.append(l[1])
        var1.append(math.sqrt(abs(l[2])) * 3.0)
        var2.append(-math.sqrt(abs(l[2])) * 3.0)

plt.plot(time, residue, 'r')
plt.plot(time, var1, 'g')
plt.plot(time, var2, 'g')
plt.title("Residue vs time plot")
plt.xlabel("Time")
plt.ylabel("Residue")
plt.show()
