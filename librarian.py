import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#def recreate_fft(data):

#def convert_to_float(data):

df1 = pd.read_csv("./fdtd_1d/saved_data/P_one_lambda.csv", skiprows=[0], sep=',', header=None)
df1 = df1.T
P = np.array(df1[0])

for elem, index in zip(P, range(0, 60000)):
    P[index] = float(elem)

fig, axes = plt.subplots()
axes.plot(P)
plt.show()


