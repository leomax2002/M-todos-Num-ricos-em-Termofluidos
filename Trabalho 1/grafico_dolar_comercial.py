import numpy as np
import matplotlib.pyplot as plt

dolar = np.array([4.752, 4.804, 4.789, 4.779, 4.796, 4.874, 4.890, 4.916, 4.987, 5.115, 5.133, 5.027, 5.027, 5.143, 5.187, 5.155, 5.177, 5.229, 5.253, 5.233])

dias = np.array ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

print(type(dolar))

graf = plt.figure()
ax = graf.add_subplot()

graf.suptitle('Dolar Comercial ao longo de 20 dias', fontsize = 12, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'Dolar', fontsize = 18, usetex = False)
ax.set_xlabel(r'Dias', fontsize = 18, usetex = False)

plt.plot(dias, dolar, '-r')
plt.plot(dias, dolar, '--b')

plt.show()

plt.savefig('dolar.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')