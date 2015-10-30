﻿#Расчёт пропеллера по импульсной теории

#Подключение модулей

import sys
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from pylab import load
from pylab import save

from scipy import interpolate
from math import pi

cosd = lambda angle: np.cos(np.deg2rad(angle))
sind = lambda angle: np.sin(np.deg2rad(angle))

#Функция чтения файла с данными профиля
def foil_data(filename, dtype, separator=';'):
	''' Считывает файл построчно с учётом разделителя ;
	Файл должен содержать колонки RE CX CY MZ, разделённые ;
	'''
	print('Source file with foil data is ', filename)
	cast = np.cast
	data = [[] for dummy in range(len(dtype))]
	for line in open(filename, 'r'):
		fields = line.strip().split(separator)
		for i, number in enumerate(fields):
			data[i].append(number)
	for i in range(len(dtype)):
		if data[i][0] == 'RE' or data[i][0] == 're' or data[i][0] == 'Re':
			print ('RE is OK')
		elif data[i][0] == 'ALPHA' or data[i][0] == 'alpha' or data[i][0] == 'Alpha':
			print ('ALPHA is OK')
		elif data[i][0] == 'CX' or data[i][0] == 'Cx' or data[i][0] == 'cx':
			print ('CX is OK')
		elif data[i][0] == 'CY' or data[i][0] == 'Cy' or data[i][0] == 'cy':
			print ('CY is OK')
		elif data[i][0] == 'MZ' or data[i][0] == 'Mz' or data[i][0] == 'mz':
			print ('MZ is OK')
		else:
			print ('ERROR in airfoil data file opening\n')
			break;
		data[i].pop(0)
	for i in range(len(dtype)):
		data[i] = cast[dtype[i]](data[i])
	return np.rec.array(data, dtype=dtype)

def prop_data(filename, dtype, separator=';'):
	''' Считывает файл построчно с учётом разделителя ;
	Файл должен содержать колонки r/R    c/R     beta, разделённые  
	'''
	print 'Source file with foil data is %s' %(filename)
	cast = np.cast
	data = [[] for dummy in range(len(dtype))]
	for line in open(filename, 'r'):
		fields = line.strip().split(separator)
		for i, number in enumerate(fields):
			data[i].append(number)
	for i in range(len(dtype)):
		if data[i][0] == 'r/R':
			print ('r/R is OK')
		elif data[i][0] == 'c/R':
			print ('c/R is OK')
		elif data[i][0] == 'beta':
			print ('beta is OK')
		else:
			print ('ERROR in geom file opening\n')
			break;
	data[0].pop(0); data[1].pop(0); data[2].pop(0)
	for i in range(len(dtype)):
		data[i] = cast[dtype[i]](data[i])
	return np.rec.array(data, dtype=dtype)

def CxCyFromFile(RE, ALPHA, foilDataArray):
	#dataFile = foil_data(sourceFile, polarDataDescr)
	if RE <= np.min(foilDataArray.RE):
		RE_MAX = np.min(foilDataArray.RE)
		RE_MIN = np.min(foilDataArray.RE)
		k_min = 0.5; k_max = 0.5
	elif RE >= np.max(foilDataArray.RE):
		RE_MAX = np.max(foilDataArray.RE)
		RE_MIN = np.max(foilDataArray.RE)
		k_min = 0.5; k_max = 0.5
	else:
		#print "Uniq elems = %s" %(np.unique(foilDataArray.RE))
		RE_SET = np.unique(foilDataArray.RE)
		RE_MAX = np.take(RE_SET, np.where(RE_SET>RE))[0][0]
		RE_MIN = np.take(RE_SET, np.where(RE_SET<RE))[0][0]	
		k_min = 1.0*(RE_MAX - RE)/(RE_MAX-RE_MIN)
		k_max = 1.0*(RE - RE_MIN)/(RE_MAX-RE_MIN)

	#print "RE_max = %s \n RE_min = %s" %(RE_MAX, RE_MIN)
	#print "k_max = %s \n k_min = %s" %(k_max, k_min)
	idx_max = np.where(foilDataArray.RE==RE_MAX)
	idx_min = np.where(foilDataArray.RE==RE_MIN)
	ALPHA_max = np.take(foilDataArray.ALPHA, idx_max)
	CX_max = np.take(foilDataArray.CX, idx_max)
	CY_max = np.take(foilDataArray.CY, idx_max)
	ALPHA_min = np.take(foilDataArray.ALPHA, idx_min)
	CX_min = np.take(foilDataArray.CX, idx_min)
	CY_min = np.take(foilDataArray.CY, idx_min)
	
	CLiftIntMax = lambda alpha: np.interp(alpha, ALPHA_max[0], CY_max[0])
	CDragIntMax = lambda alpha: np.interp(alpha, ALPHA_max[0], CX_max[0])
	CLiftIntMin = lambda alpha: np.interp(alpha, ALPHA_min[0], CY_min[0])
	CDragIntMin = lambda alpha: np.interp(alpha, ALPHA_min[0], CX_min[0])
	
	return [CDragIntMax(ALPHA)*k_max + CDragIntMin(ALPHA)*k_min,
	CLiftIntMax(ALPHA)*k_max + CLiftIntMin(ALPHA)*k_min]
	
#Считываем Данные поляры

if len(sys.argv)> 0:
	polarDataDescr = np.dtype([('RE', 'int'),('ALPHA', 'float'),  ('CX', 'float'), ('CY', 'float'), ('MZ', 'float')])
	polar_file = sys.argv[1]
	foilData = foil_data(polar_file, polarDataDescr)
	CLift = lambda re, alpha: CxCyFromFile(re, alpha, foilData)[1]
	CDrag = lambda re, alpha: CxCyFromFile(re, alpha, foilData)[0]
else:
	print('No foil data file has been specified ')
	sys.exit()

Ntab = 9
N = 3008
W = np.pi/30*N
V0 = np.array([0.01, 2, 4, 5, 6, 6.5, 7, 7.5, 8])
J = V0/(N/60*10*0.025)
print "V0 = %s" %V0
print "J = %s" %J
print "W*R = %s" %(W*10*0.025*0.5)
angle = np.zeros(Ntab)
Cx_interp = np.zeros(Ntab)

#Парметры винта
if len(sys.argv)> 1:
	geom_file = sys.argv[2]
	geomDataDescr = np.dtype([('r_', 'float'),('c_', 'float'),  ('beta', 'float')])
	geomData = prop_data(geom_file, geomDataDescr )
	#print geomData.r_
	b = lambda r_i: np.interp(r_i, geomData.r_, geomData.c_) #Функция вычисления хорды на основе файла геометрии
	sigma = lambda r_i: K_l*b(r_i)/(np.pi)
	fPhi = lambda r_i: np.interp(r_i, geomData.r_, geomData.beta) #функция для расчёта угла установки от относительного радиуса 
else:
	print('No propeller geometry source file has been specified ')
	sys.exit()
	
#----------------------------------------------#			
D = 10*0.0254	#дюйм*коэфф = м,   диаметр винта
R_l = 0.5*D		#м	
H = 7*0.0254	#дюйм*коэфф = м,   шаг винт
b_07 = 0.214*R_l	#м	хорда на 0,7
eta_l = 1	# 	сужение лопасти в плане;
K_l = 2		#	количество лопастей НВ;
r_0 = 0.15	#	относительный радиус втулки, с которого начинается рабочее сечение лопасти;
Nrad = 11	#	количество сечений по r
a_inf = 5.6
a = 341.4
nu = 1.5e-5
#----------------------------------------------#			

#функция для расчёта коэффициентa заполнения от относительной хорды 	

r = np.zeros(Nrad)

#Вычисляем отн. радиус, угол установки, коэффициент заполнения 
for i in range(1, Nrad):
	r[i]=r[i-1]+1.0/(Nrad-1.0)

Mach = lambda omega, R, soundVel: omega*R/soundVel
Reynolds = lambda omega, R, chord, visc: omega*R*chord/visc

Mch = np.zeros((Nrad, Ntab))
Re = np.zeros((Nrad, Ntab))

#Вычисляем М и Re по оборотам и по размаху
for i in range(0, Nrad):
	for j in range(0,Ntab):
		Mch[i][j] = Mach( W, R_l*r[i], a)
		Re[i][j] = Reynolds( W, r[i]*R_l, b(r[i])*R_l, nu)

print Re[1,0]
V_1 = np.zeros((Nrad, Ntab))
beta_e = np.zeros((Nrad, Ntab))
alpha_e = np.zeros((Nrad, Ntab))

#Функции для интерполяции происходящего
CLiftTest = lambda alpha: np.interp(alpha, ALPHA_work[0], CY_work[0])
CDragTest = lambda alpha: np.interp(alpha, ALPHA_work[0], CX_work[0])

print "CLIFT LAMBDA interpolation test"
print "RE = %s, Cy = %s" %(75000, CLift(75000, 5.0))
print "RE = %s, Cx = %s" %(75000, CDrag(75000, 5.0))

def solve_V1(CyFunct, CxFunct, Re_e, V_0, sigmaVar, rVar , phi_e):
	if rVar == 0:
		return [0,0]
	beta_e = 0
	V_1_old = 0
	for i in range(0,100):
		alpha_e = phi_e - beta_e
		Cy_beta = CyFunct(Re_e, alpha_e)*cosd(beta_e) - CxFunct(Re_e, alpha_e)*sind(beta_e)
		C = 0.125*Cy_beta*sigmaVar/rVar
		V_1 = (V_0 + np.sqrt( V_0**2 + 4*(1-C)*C*rVar**2))*0.5/(1-C)
		beta_new = np.rad2deg(np.arctan(V_1/rVar))
		if np.abs(beta_new-beta_e)/(beta_new)<0.001:
			print "V_1 converged for n=%i iterations" %i
			print "Delta V_1 = %s " %(np.abs(V_1_old - V_1)/V_1_old)
			break
		else:
			beta_e = beta_new
			V_1_old = V_1
	return [V_1, alpha_e]

V_0 = 0.0

#Расчёт скоростей набегания
for i in range(0, Nrad):
	for j in range(0,Ntab):
		if W != 0:
			[V_1[i,j], alpha_e[i,j]] = solve_V1(CLift, CDrag, Re[i,j], V0[j]/(R_l*W), sigma(r[i]), r[i] , fPhi(r[i]))
			beta_e[i,j] = fPhi(r[i])-alpha_e[i,j]
		else:
			[V_1[i,j], alpha_e[i,j] ] = [0, 0]
			beta_e[i,j] = 0

plt.figure(1)
plt.subplot(411)
plt.plot(fPhi(r))
#plt.axis([0, 9000, 0, 1000])
plt.ylabel('phi_e, grad')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(412)
graphs = plt.plot(r, alpha_e[:,:])
#plt.axis([0, 9000, 0, 1000])
plt.ylabel('alpha_e, grad')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(413)
plt.plot(r, beta_e[:,:])
#plt.axis([0, 9000, 0, 10])
plt.ylabel('beta_e, grad')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(414)
plt.plot(r, V_1[:,:])
#plt.axis([0, 9000, 0, 5])
plt.xlabel('Relative radius')
plt.ylabel('Inductive velocity')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

#plt.show()
plt.savefig('angles.png', bbox_inches=0)

dCt_dr = np.zeros((Nrad, Ntab))
dCt_dr_v = np.zeros((Nrad, Ntab))

for i in range(1, Nrad):
	for j in range(0,Ntab):
		dCt_dr[i,j] = CLift(Re[i,j], alpha_e[i,j])*sigma(r[i])*0.25*(r[i]+r[i-1])**2*(r[i]-r[i-1])
		dCt_dr_v[i,j] = V_1[i,j]*(V_1[i,j]-V0[j]/(R_l*W))*(r[i]+r[i-1])*0.5*(r[i]-r[i-1])

Ct_c = np.sum(dCt_dr, axis = 0)*np.pi**3/8
Ct_c_v = np.sum(dCt_dr_v, axis = 0)*np.pi**3
#Концевые потери
B = 1 - 4*Ct_c/K_l
Ct = B*Ct_c

print('Ct* list\n', Ct_c)	
print('Ct_v* list\n', Ct_c_v)	

print('B list\n', B)	
print('Ct list\n', Ct)

dmi_dr = np.zeros((Nrad, Ntab))
dmr_dr = np.zeros((Nrad, Ntab))
dmi_dr_v = np.zeros((Nrad, Ntab))
dmr_dr_v = np.zeros((Nrad, Ntab))

for i in range(1, Nrad):
	for j in range(0,Ntab):
		dmi_dr[i,j] = dCt_dr[i,j]*V_1[i,j]*(r[i]-r[i-1])
		dmr_dr[i,j] = CDrag(Re[i,j], alpha_e[i,j])*sigma(r[i])*0.125*(r[i-1]+r[i])**3*(r[i]-r[i-1])
		
		dmi_dr_v[i,j] = dCt_dr_v[i,j]*V_1[i,j]*(r[i]-r[i-1])
		Cx_beta = CLift(Re[i,j], alpha_e[i,j])*sind(beta_e[i,j]) + CDrag(Re[i,j], alpha_e[i,j])*cosd(beta_e[i,j])
		dmr_dr_v[i,j] = Cx_beta*sigma(r[i])*(r[i]**2 + V_1[i,j]**2)*(r[i]+r[i-1])*0.5*(r[i]-r[i-1])

mi = np.sum(dmi_dr, axis = 0)
mr = np.sum(dmr_dr, axis = 0)
mk = (mi + mr)*np.pi**4/8
print "Mk"
print mk
print('dmi/dr \n', mi)	
print('dmr/dr \n', mr)	
#print('dmr/dr \n', dmr_dr)	
mi_v = np.sum(dmi_dr_v, axis = 0)
mr_v = np.sum(dmr_dr_v, axis = 0)
mk_v = (mi_v + mr_v)*np.pi**4/8
print "Mk_v"
print mk_v
print('dmi/dr \n', mi_v)	
print('dmr/dr \n', mr_v)	

nu0 = Ct_c_v/mk_v*J
print('Nu0 list\n', nu0)	

plt.figure(2)
plt.subplot(311)
graphs = plt.plot(J, Ct_c_v)
plt.ylabel('Ct')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(312)
plt.plot(J, mk_v)
plt.ylabel('Cp')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(313)
plt.plot(J, nu0)
plt.xlabel('Advance ratio')
#plt.axis([0, 1, 0, 1])
plt.ylabel('Efficiency')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.show()
plt.savefig('results.png', bbox_inches=0)

#Сохраняем данные в файл
outfile=open('N'+str(N)+'.dat', 'w')
outfile.write('J \t Cp \t Ct \t eta\n')

for j, cp, ct, eta in zip(J, Ct_c_v, mk_v, nu0):
	outfile.write(str(j) + '\t'+str(ct) + '\t' + str(cp) + '\t' + str(eta) + '\n')
outfile.close()

#Находим мощность на уровне моря

rho = 1.226
'''
Nv = 0.5*mk*rho*(w*R_l)**3*np.pi*R_l**2
G = 0.5*Ct*rho*(w*R_l)**2*np.pi*R_l**2
m = G/9.81
Mn = Nv/(w*R_l)
print('Nv list\n', Nv)	
print('G list\n', G)	
print('m list\n', m)	
print('m list\n', Mn)	

plt.figure(1)
plt.subplot(311)
plt.plot(n, N, 'b',n, Nv, 'r')
plt.axis([0, 9000, 0, 1000])
plt.xlabel('n, rev/min')
plt.ylabel('N, Nv, W')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(312)
plt.plot(n, Mn, 'r', n , M, 'b')
plt.axis([0, 9000, 0, 10])
plt.xlabel('n')
plt.ylabel('M, Mn, N/m')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

plt.subplot(313)
plt.plot(n, m, 'b')
plt.axis([0, 9000, 0, 5])
plt.xlabel('n')
plt.ylabel('G, kg')
plt.grid(True)
plt.axhline(0, color='black', lw=1)

#plt.show()
plt.savefig('moments.png', bbox_inches=0)
'''