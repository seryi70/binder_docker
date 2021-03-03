# coding: utf-8
"""
Created on SunMay 15 18:23:58 2016
@author: Spencer Pringle @Rochester Institute of Technology
"""

from pylab import *
from PyQt4 import QtGui , QtCore
from matplotlib.backends.backend_qt4agg import (FigureCanvasQTAgg as FigureCanvas ,NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import sys
import GUI_NEW

class ExampleApp(QtGui.QMainWindow , GUI_NEW.Ui_MainWindow):
	def __init__(self):
		super(self.__class__, self).__init__()
		self.setupUi(self)
		self.ModelingProgress.setValue(0)
		self.Metal1Box.addItems(['SrRuO_3 ','Co','Al','Metal 1'])
		self.Metal1Box.setCurrentIndex(0)
		self.Metal2Box.addItems(['PZT','La_(0.67)Sr_(0.33)MnO_3 ','p+ Silicon', 'Metal 2'])
		self.Metal2Box.setCurrentIndex(0)
		self.FerroBox.addItems(['BaTiO_3 ','Al-HfO_2 ', 'Ferro'])
		self.FerroBox.setCurrentIndex(0)
		self.RunModeling.clicked.connect(self.Model)
		self.Metal1Box.currentIndexChanged.connect(self.Metal1Update)
		self.Metal2Box.currentIndexChanged.connect(self.Metal2Update)
		self.FerroBox.currentIndexChanged.connect(self.FerroelectricUpdate)
		self.Metal1Box.setCurrentIndex(2)
		self.Metal2Box.setCurrentIndex(2)
		self.FerroBox.setCurrentIndex(1)

	def Model(self):
		self.ModelingProgress.setValue(0)
		self.FrontOutput.clear()
		self.PotentialOutput.clear()
		global sigma_p , sigma_s , q, epsilon_0 , phi, sigma , k_1 , k_2 , a_0 , x,y, meshSpace , potentialMesh , phi_2 , phi_1 , positiveMesh ,negativeMesh
		global e0, m0, kT, kb, h, m_0
		global d, E_a , epsilon_f , P, h, h_eV
		global delta_1 , a_1 , E_f_1
		global delta_2 , a_2 , E_f_2
		
		
		global G_2 , G_1 , x_1 , x_2 , phi_bar_pos , phi_bar_neg , x_1_index ,x_2_index , m_0 , A_tun , J_0 , J, v_app , phi_diff
		global Ath_pos , Jth_pos , Ath_neg , Jth_neg , phi_prime , phi_prime_pos ,phi_prime_neg
		global Atun_pos , J0_pos , Jtun_pos , Atun_neg , J0_neg , Jtun_neg , V, T,R_pos , R_neg , m_star
		global Ran, numchild , i, w, items , WellDir , WellFermi , Bar_Length_Pos, Bar_Length_Neg , x_1_pos , x_2_pos , V, h_bar
		
		epsilon_0 =8.854E-12#free space perm [F/m]
		T=300 #Temperature (K)
		#Non-adjustable Constants
		e0=8.854E-14 #Permittivity of free space [F/cm]
		m_0 =9.11E-31 #Electron rest mass (kg)
		h=6.626E-34 #Planck's constant (m^2*kg/sec)
		q=1.6E-19 #Electron charge (C)
		kb=1.38E-23 #Boltzman constant(J/K)
		kT=(kb/q)*T #[eV]
		a_0 =5.291E-11
		h=6.626E-34
		h_bar=h/(2*pi)
		h_eV =4.13E-15
		meshSpace =2/1000
		v_app =.1
		
		#Initialize Ferroelectric Material Values
		epsilon_f=self.FerroDielectricConst.value()
		
		
		E_a =self.FerroEA.value()
		m_star=self.FerroEMass.value()
		P=self.FerroPolarization.value()
		d=self.FerroThick.value()*1E-9
		
		#Initialize Metal 1 Values
		E_f_1 =self.Metal1FermiEnergy.value()
		if self.Metal1ScreenLength.value()==0:
			delta_1 =( E_f_1 /(4*pi*q*self.Metal1Lattice.value()*1E-10))**(1/2)
		else:
			delta_1=self.Metal1ScreenLength.value()*1E-9
		
		#Initialize Metal 2 Values
		#if self .Metal2Box.currentIndex()==2:
		E_f_2 =self.Metal2FermiEnergy.value()
		#else:
		# E_f_2=self .Metal2FermiEnergy.value()
		if self.Metal2ScreenLength.value()==0:
			delta_2 =( E_f_2 /(4*pi*q*self.Metal2Lattice.value()*1E-10))**(1/2)
		else:
			delta_2=self.Metal2ScreenLength.value()*1E-9
		
		#Calculate and store charges and wave vector magnitudes...
		sigma_p=P*1E-2
		sigma_s=sigma_p*d/( epsilon_f*(delta_1+delta_2)+d)
		V=np.linspace(-1,1,201)
		self.Potential()
		x_1_pos = x_1
		x_2_pos = x_2
		positiveMesh=potentialMesh
		phi_prime_pos=phi_prime
		Ath_pos=(4*pi*m_0*(kb**2)*q)/((h**3))
		Jth_pos = Ath_pos*1E-4*(T**2)*exp(-phi_prime_pos/kT)*(1-exp(-abs(V)/kT))
		phi_diff=phi_1-phi_2
		Bar_Length_Pos=x_2-x_1
		phi_bar_pos =mean(y[ x_1_index : x_2_index ])
		Atun_pos=((2*m_0*q*m_star)**(1/2))*(Bar_Length_Pos)*(2e-9)/h_bar
		#Atun_pos=(4*pi*(x_2-x_1)*1e-9*(2*m0)**(1/2))/h
		#J0_pos=(q)/(2*pi*h*(((x_2-x_1)*1E-9)**2))
		J0_pos =(6.08e8)/(( Bar_Length_Pos)**2)
		Jtun_pos=J0_pos*((phi_bar_pos )*exp(-Atun_pos*((phi_bar_pos )**(1/2)))-(phi_bar_pos +V)*exp(-Atun_pos*((phi_bar_pos +V)**(1/2))))
		R_pos=np.zeros(size(V))
		for i in range(0, size(V)):
			if Jtun_pos[i]!=0 and Jtun_pos[i]!='nan':
				R_pos[i]=V[i]/ Jtun_pos[i]
			else:
				R_pos[i]='nan'
		np.nanmean(abs(R_pos))
		Jtun_pos=abs(Jtun_pos)
		Itun_pos=Jtun_pos[np.where(V==.5)[0]]*(250e-7)**2
		Itun_pos_1_250 =Jtun_pos [[110][0]]*(250e-7)**2
		Itun_pos_2_250 =Jtun_pos [[120][0]]*(250e-7)**2
		Itun_pos_1_500 =Jtun_pos [[110][0]]*(500e-7)**2
		Itun_pos_2_500 =Jtun_pos [[120][0]]*(500e-7)**2
		Rtun_pos_250 =.5/ Itun_pos
		Rtun_pos_1_250 =.1/ Itun_pos_1_250
		Rtun_pos_1_500 =.1/ Itun_pos_1_500
		Rtun_pos_2_250 =.2/ Itun_pos_2_250
		Rtun_pos_2_500 =.2/ Itun_pos_2_500
		rho_pos =.5/ Jtun_pos[np.where(V==.5)[0]]
		rho_pos_1 =.1/ Jtun_pos [[110][0]]
		self.FrontOutput.insertPlainText('Under Positive Polarization (towardMetal 1), at .5V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_pos[np.where(V==.5)[0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_pos[np.where(V==.5)[0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_pos +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Under Positive Polarization (towardMetal 1), at .2V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_pos [[120][0]]+' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_pos[[120][0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_pos_2_250 +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Under Positive Polarization (towardMetal 1), at .1V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_pos [[110][0]]+' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_pos[[110][0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_pos_1_250 +' A')
		self.FrontOutput.insertPlainText('\n')
		
		P=-P
		sigma_p=P*1E-2
		sigma_s=sigma_p*d/( epsilon_f*(delta_1+delta_2)+d)
		
		
		self.Potential()
		x_1_neg = x_1
		x_2_neg = x_2
		negativeMesh=potentialMesh
		phi_prime_neg=phi_prime
		Ath_neg=(4*pi*m_0*(kb**2)*q)/((h**3))
		Jth_neg = Ath_neg*1E-4*(T**2)*exp(-phi_prime_neg/kT)*(1-exp(-abs(V)/kT))
		
		Bar_Length_Neg=x_2-x_1
		phi_bar_neg =mean(y[ x_1_index : x_2_index ])
		Atun_neg=((2*m_0*q*m_star)**(1/2))*(Bar_Length_Neg)*(2e-9)/h_bar
		J0_neg =(6.08e8)/(( Bar_Length_Neg)**2)
		Jtun_neg=J0_neg*((phi_bar_neg )*exp(-Atun_neg*((phi_bar_neg )**(1/2)))-((phi_bar_neg )+V)*exp(-Atun_neg*(((phi_bar_neg )+V)**(1/2))))
		R_neg=np.zeros(size(V))
		for i in range(0, size(V)):
			if Jtun_neg[i]!=0 and Jtun_neg[i]!='nan':
				R_neg[i]=V[i]/ Jtun_neg[i]
			else:
				R_neg[i]='nan'
		Jtun_neg=abs(Jtun_neg)
		Itun_neg=Jtun_neg[np.where(V==.5)[0]]*(250e-7)**2
		Itun_neg_1_250 =Jtun_neg [[110][0]]*(250e-7)**2
		Itun_neg_1_500 =Jtun_neg [[110][0]]*(500e-7)**2
		Itun_neg_2_250 =Jtun_neg [[120][0]]*(250e-7)**2
		Itun_neg_2_500 =Jtun_neg [[120][0]]*(500e-7)**2
		
		
		Rtun_neg_250 =.5/ Itun_neg
		Rtun_neg_1_250 =.1/ Itun_neg_1_250
		Rtun_neg_1_500 =.1/ Itun_neg_1_500
		Rtun_neg_2_250 =.2/ Itun_neg_2_250
		Rtun_neg_2_500 =.2/ Itun_neg_2_500
		rho_neg =.5/ Jtun_neg[np.where(V==.5)[0]]
		rho_neg_1 =.1/ Jtun_neg [[110][0]]
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Under Negative Polarization (towardMetal 2), at .5V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_neg[np.where(V==.5)[0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_neg[np.where(V==.5)[0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_neg +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Under Negative Polarization (towardMetal 2), at .2V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_neg [[120][0]]+' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_neg[[120][0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_neg_2_250 +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (500nm x 500nm)='+'%.2e'%Itun_neg_2_500 +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Under Negative Polarization (towardMetal 2), at .1V bias:')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J th='+'%.2e'%Jth_neg [[110][0]]+' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t J tun='+'%.2e'%Jtun_neg[[110][0]] +' A/cm^2')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (250nm x 250nm)='+'%.2e'%Itun_neg_1_250 +' A')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\t I tun (500nm x 500nm)='+'%.2e'%Itun_neg_1_500 +' A')
		self.FrontOutput.insertPlainText('\n')
		
		phi_diff=abs(phi_bar_pos-phi_bar_neg )
		if Jtun_neg[np.where(V==.5)[0]]>Jtun_pos[np.where(V==.5)[0]]:
			Jtun_Ratio=Jtun_neg[np.where(V==.5)[0]]/ Jtun_pos[np.where(V==.5)[0]]
			LRSR=np.nanmean(abs(R_neg))
			HRSR=np.nanmean(abs(R_pos))
		else:
			Jtun_Ratio=Jtun_pos[np.where(V==.5)[0]]/ Jtun_neg[np.where(V==.5)[0]]
			LRSR=np.nanmean(abs(R_pos))
			HRSR=np.nanmean(abs(R_neg))
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Difference in Average PotentialBarriers:\n')
		self.FrontOutput.insertPlainText('\t'+'%.3f'% phi_bar_pos +'-'+'%.3f'%phi_bar_neg +'='+'%.3f'%phi_diff + ' eV\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('Ratio of Tunnel Current at .5V Biasin HRS vs. LRS:\n')
		self.FrontOutput.insertPlainText('\t'+'%.3f'%Jtun_Ratio+'\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistances (250nm x250nm, at .5V):\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% Rtun_pos_250 +' Ohm\n\t Neg:'+'%.3e'% Rtun_neg_250 +' Ohm\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistivity @ .5V:\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'%rho_pos +' Ohm cm^2\n\t Neg:'+'%.3e'%rho_neg +' Ohm cm^2\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistivity @ .1V:\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% rho_pos_1 +' Ohm cm^2\n\t Neg:'+'%.3e'% rho_neg_1 +' Ohm cm^2\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistances (250nm x250nm, at .2V):\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% Rtun_pos_2_250 +'Ohm\n\t Neg:'+'%.3e'% Rtun_neg_2_250 +' Ohm\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistances (500nm x500nm, at .2V):\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% Rtun_pos_2_500 +'Ohm\n\t Neg:'+'%.3e'% Rtun_neg_2_500 +' Ohm\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistances (250nm x250nm, at .1V):\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% Rtun_pos_1_250 +'Ohm\n\t Neg:'+'%.3e'% Rtun_neg_1_250 +' Ohm\n')
		self.FrontOutput.insertPlainText('\n')
		self.FrontOutput.insertPlainText('LRS and HRS Resistances (500nm x500nm, at .1V):\n')
		self.FrontOutput.insertPlainText('\t Pos:'+'%.3e'% Rtun_pos_1_500 +'Ohm\n\t Neg:'+'%.3e'% Rtun_neg_1_500 +' Ohm\n')
		self.FrontOutput.insertPlainText('\n')
		Jtunratiomax = max(np.nan_to_num (np.divide(Jtun_pos , Jtun_neg)))
		self.FrontOutput.insertPlainText('Max Current Ratio: '+'%.3e'%Jtunratiomax)
		self.closeall()
		fig1 = Figure()
		ax1f1 = fig1.add_subplot(111)
		ax1f1.plot(positiveMesh[0,:], positiveMesh[1,:])
		ax1f1.plot(negativeMesh[0,:], negativeMesh[1,:])
		ExampleApp.addmpl(self,fig1)
		fig2=Figure()
		ax1f2 = fig2.add_subplot(111)
		ax1f2.plot(V, Jth_pos)
		ax1f2.plot(V, Jth_neg)
		ExampleApp.addmpl(self,fig2)
		fig3=Figure()
		ax1f3 = fig3.add_subplot(111)
		ax1f3.plot(V, Jtun_pos)
		ax1f3.plot(V, Jtun_neg)
		ax1f3.set_yscale('log')
		ExampleApp.addmpl(self,fig3)
		self.PotentialOutput.insertPlainText('V\tJtunpos\tJtunneg\n')
		for i in range(0, size(V)):
			self.PotentialOutput.insertPlainText('%.3f'%V[i] + '\t' + '%.3e'%Jtun_pos[i] + '\t' + '%.3e'%Jtun_neg[i] + '\n')
		self.PotentialOutput.insertPlainText('\n')
		self.PotentialOutput.insertPlainText('X\tV1pos\tV2neg\n')
		for p in range(0, size(x)):
			self.PotentialOutput.insertPlainText('%.3f'%x[p] + '\t' + '%.3f'%positiveMesh[1,p] + '\t' + '%.3f'%negativeMesh[1,p]+ '\n')

	def unfill(self):
		global numchild , i, w, items
		numchild=self.PotentialMPlay.count()
		items = (self.PotentialMPlay.itemAt(i) for i in range(self.PotentialMPlay.count()))
		for w in items:
			try:
				self.PotentialMPlay.removeWidget(w.widget())
			except AttributeError:
				pass
		try:
			self.PotentialMPlay.removeItem(w.item())
		except AttributeError:
			pass

	def closeall(self): #solution found on http://python.6.x6.nabble.com/Completely-removing-items-from-a-layout-td1924202.html
		def deleteItems(layout):
			if layout is not None:
				while layout.count():
					item = layout.takeAt(0)
					widget = item.widget()
					if widget is not None:
						widget.deleteLater()
					else:
						deleteItems(item.layout())
		deleteItems(self.PotentialMPlay)

	def Potential(self):
		global P, x, index, y, x_prime, potentialMesh, delta_1, delta_2, E_f_1, E_f_2, phi_2, phi_1, x_1,x_2, x_1_index, x_2_index,phi_prime, WellDir, WellFermi
		
		x = np.linspace(-5,5+(d*1E9),(5+(d*1E9))/meshSpace)
		y = np.zeros(size(x))
		index = 0
		flag = 0
		y_flag =0
		
		if delta_1<delta_2:
			WellFermi= E_f_2
			WellDir=-1
		
		
		else:
			WellFermi= E_f_1
			WellDir=1
		
		for x_prime in x:
			if x_prime<=0:
				y[index]=( sigma_s*delta_1/epsilon_0)*exp(-abs(x_prime*1E-9)/( delta_1))
			elif x_prime>=(d*1E9):
				y[index]=(-sigma_s*delta_2/epsilon_0)*exp(-abs((x_prime*1E-9)-d)/( delta_1))
			if x_prime<=0:
				phi_1=y[index]
			if x_prime>=(d*1E9) and flag==0:
				phi_2=y[index]
				flag=1
			index=index+1
			self.ModelingProgress.setValue(index*50/size(x))
		index = 0
		for x_u in x:
			if x_u>0 and x_u<(d*(1E9)):
				y[index]=phi_1-((x_u /(d*1E9))*(phi_1-phi_2))
			index=index+1
			self.ModelingProgress.setValue(50+(index*25/size(x)))
		index = 0
		for x_u in x:
			if x_u>0 and x_u<(d*(1E9)):
				y[index]=y[index]+( E_f_1-E_a )+(( x_u /(d*1E9))*(E_f_2-E_f_1 ))
			if x_u<=0:
				phi_1=y[index]
			if x_u>=(d*1E9) and flag==0:
				phi_2=y[index]
				flag=1
			index=index+1
			self.ModelingProgress.setValue(50+(index*25/size(x)))
		phi_prime=max(y)
		index=0
		for x_u in x:
			if y_flag ==0 and y[index]>phi_prime /10:
				x_1 = x_u
				x_1_index =index
				y_flag =1
			if y_flag ==1 and y[index]<phi_prime /10:
				x_2 = x_u
				x_2_index =index-1
				y_flag =2
			index=index+1
			self.ModelingProgress.setValue(75+(index*25/size(x)))
		potentialMesh=vstack((x,y))

	def addmpl(self, fig):
		for clos in range(self.PotentialMPlay.count()):
			try:
				self.canvas.close()
				self.toolbar.close()
			except AttributeError:
				pass
		self.canvas = FigureCanvas(fig)
		self.PotentialMPlay.addWidget(self.canvas)
		self.canvas.draw()
		self.toolbar = NavigationToolbar(self.canvas ,self.verticalWidget_6 , coordinates=True)
		self.PotentialMPlay.addWidget(self.toolbar)

	def Metal1Update(self):
		if self.Metal1Box.currentIndex()==0:
			self.Metal1ScreenLength.setValue(.6)
			self.Metal1FermiEnergy.setValue(1.5)
			self.Metal1Lattice.setValue(0)
		elif self.Metal1Box.currentIndex()==1:
			self.Metal1ScreenLength.setValue(.05)
			self.Metal1FermiEnergy.setValue(5)
			self.Metal1Lattice.setValue(0)
		elif self.Metal1Box.currentIndex()==2:
			self.Metal1ScreenLength.setValue(.06)
			self.Metal1FermiEnergy.setValue(4.08)
			self.Metal1Lattice.setValue(0)

	def FerroelectricUpdate(self):
		if self.FerroBox.currentIndex()==0:
			self.FerroDielectricConst.setValue(2000)
			self.FerroEA.setValue(2.5)
			self.FerroPolarization.setValue(20)	
			self.FerroThick.setValue(2)
			self.FerroEMass.setValue(1)
		elif self.FerroBox.currentIndex()==1:
			self.FerroDielectricConst.setValue(40)
			self.FerroEA.setValue(2.0)#http://e-citations.ethbib.ethz.ch/view/pub:28311
			self.FerroPolarization.setValue(10)
			self.FerroThick.setValue(2)#http://e-citations.ethbib.ethz.ch/view/pub:28311
			self.FerroEMass.setValue(.11)

	def Metal2Update(self):
		if self.Metal2Box.currentIndex()==0:
			self.Metal2ScreenLength.setValue(.07)
			self.Metal2FermiEnergy.setValue(3.5)
			self.Metal2Lattice.setValue(0)
		elif self.Metal2Box.currentIndex()==1:
			self.Metal2ScreenLength.setValue(.1)
			self.Metal2FermiEnergy.setValue(4.8)#FromAbuwasib 15
			self.Metal2Lattice.setValue(0)
		elif self.Metal2Box.currentIndex()==2:
			self.Metal2ScreenLength.setValue(.4)
			self.Metal2FermiEnergy.setValue(4.85)
			self.Metal2Lattice.setValue(0)

def main():
	if QtCore.QCoreApplication.instance() != None:
		app = QtCore.QCoreApplication.instance()
	else:
		app = QtGui.QApplication(sys.argv) #Anew instance of QApplication
		form = ExampleApp() #We set the form to be ourExampleApp (design)
		form.show() #Show the form
		app.exec_() #and execute the app
if __name__ == '__main__': # if we're running file directly and not importing it
	main() #run the main function

try:
	_fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
	def _fromUtf8(s):
		return s

try:
	_encoding = QtGui.QApplication.UnicodeUTF8
	def _translate(context , text , disambig):
		return QtGui.QApplication.translate(context , text , disambig ,_encoding)
except AttributeError:
	def _translate(context , text , disambig):
		return QtGui.QApplication.translate(context , text , disambig)

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName( fromUtf8("MainWindow"))
		MainWindow.resize(941, 739)
		self.centralwidget = QtGui.QWidget(MainWindow)
		self.centralwidget.setObjectName( fromUtf8("centralwidget"))
		self.gridLayout = QtGui.QGridLayout(self.centralwidget)
		self.gridLayout.setObjectName( fromUtf8("gridLayout"))
		self.splitter_5 = QtGui.QSplitter(self.centralwidget)
		self.splitter_5.setOrientation(QtCore.Qt.Vertical)
		self.splitter_5.setObjectName( fromUtf8("splitter 5"))
		self.splitter = QtGui.QSplitter(self.splitter_5)
		self.splitter.setOrientation(QtCore.Qt.Horizontal)
		self.splitter.setObjectName( fromUtf8("splitter"))
		self.layoutWidget = QtGui.QWidget(self.splitter)
		self.layoutWidget.setObjectName( fromUtf8("layoutWidget"))
		self.verticalLayout_7 = QtGui.QVBoxLayout(self.layoutWidget)
		self.verticalLayout_7.setObjectName( fromUtf8("verticalLayout_7"))
		self.verticalLayout_4 = QtGui.QVBoxLayout()
		self.verticalLayout_4.setObjectName( fromUtf8("verticalLayout_4"))
		self.verticalLayout = QtGui.QVBoxLayout()
		self.verticalLayout.setObjectName( fromUtf8("verticalLayout"))
		self.groupBox = QtGui.QGroupBox(self.layoutWidget)
		self.groupBox.setObjectName( fromUtf8("groupBox"))
		self.gridLayout_3 = QtGui.QGridLayout(self.groupBox)
		self.gridLayout_3.setObjectName( fromUtf8("gridLayout_3"))
		self.verticalLayout_3 = QtGui.QVBoxLayout()
		self.verticalLayout_3.setObjectName( fromUtf8("verticalLayout_3"))
		self.horizontalLayout_2 = QtGui.QHBoxLayout()
		self.horizontalLayout_2.setObjectName( fromUtf8("horizontalLayout_2"))
		self.label = QtGui.QLabel(self.groupBox)
		self.label.setObjectName( fromUtf8("label"))
		self.horizontalLayout_2.addWidget(self.label)
		self.label_2 = QtGui.QLabel(self.groupBox)
		self.label_2.setObjectName( fromUtf8("label 2"))
		self.horizontalLayout_2.addWidget(self.label_2)
		self.label_3 = QtGui.QLabel(self.groupBox)
		self.label_3.setObjectName( fromUtf8("label 3"))
		self.horizontalLayout_2.addWidget(self.label_3)
		self.verticalLayout_3.addLayout(self.horizontalLayout_2)
		self.horizontalLayout = QtGui.QHBoxLayout()
		self.horizontalLayout.setObjectName( fromUtf8("horizontalLayout"))
		self.Metal1Box = QtGui.QComboBox(self.groupBox)
		self.Metal1Box.setObjectName( fromUtf8("Metal1Box"))
		self.horizontalLayout.addWidget(self.Metal1Box)
		self.FerroBox = QtGui.QComboBox(self.groupBox)
		self.FerroBox.setObjectName( fromUtf8("FerroBox"))
		self.horizontalLayout.addWidget(self.FerroBox)
		self.Metal2Box = QtGui.QComboBox(self.groupBox)
		self.Metal2Box.setObjectName( fromUtf8("Metal2Box"))
		self.horizontalLayout.addWidget(self.Metal2Box)
		self.verticalLayout_3.addLayout(self.horizontalLayout)
		self.gridLayout_3.addLayout(self.verticalLayout_3 , 0, 0, 1, 1)
		self.verticalLayout.addWidget(self.groupBox)
		self.verticalLayout_4.addLayout(self.verticalLayout)
		self.verticalLayout_7.addLayout(self.verticalLayout_4)
		self.splitter_4 = QtGui.QSplitter(self.layoutWidget)
		self.splitter_4.setOrientation(QtCore.Qt.Vertical)
		self.splitter_4.setObjectName( fromUtf8("splitter 4"))
		self.layoutWidget1 = QtGui.QWidget(self.splitter_4)
		self.layoutWidget1.setObjectName( fromUtf8("layoutWidget1"))
		self.verticalLayout_12 = QtGui.QVBoxLayout(self.layoutWidget1)
		self.verticalLayout_12.setObjectName( fromUtf8("verticalLayout_12"))
		self.label_14 = QtGui.QLabel(self.layoutWidget1)
		self.label_14.setObjectName( fromUtf8("label 14"))
		self.verticalLayout_12.addWidget(self.label_14)
		self.horizontalLayout_4 = QtGui.QHBoxLayout()
		self.horizontalLayout_4.setObjectName( fromUtf8("horizontalLayout_4"))
		self.Metal1Group = QtGui.QGroupBox(self.layoutWidget1)
		self.Metal1Group.setObjectName( fromUtf8("Metal1Group"))
		self.horizontalLayout_3 = QtGui.QHBoxLayout(self.Metal1Group)
		self.horizontalLayout_3.setObjectName( fromUtf8("horizontalLayout_3"))
		self.verticalLayout_9 = QtGui.QVBoxLayout()
		self.verticalLayout_9.setObjectName( fromUtf8("verticalLayout_9"))
		self.label_4 = QtGui.QLabel(self.Metal1Group)
		self.label_4.setObjectName( fromUtf8("label 4"))
		self.verticalLayout_9.addWidget(self.label_4)
		self.Metal1ScreenLength = QtGui.QDoubleSpinBox(self.Metal1Group)
		self.Metal1ScreenLength.setDecimals(5)
		self.Metal1ScreenLength.setProperty("value", 0.6)
		self.Metal1ScreenLength.setObjectName( fromUtf8("Metal1ScreenLength"))
		self.verticalLayout_9.addWidget(self.Metal1ScreenLength)
		self.label_5 = QtGui.QLabel(self.Metal1Group)
		self.label_5.setObjectName( fromUtf8("label 5"))
		self.verticalLayout_9.addWidget(self.label_5)
		self.Metal1FermiEnergy = QtGui.QDoubleSpinBox(self.Metal1Group)
		self.Metal1FermiEnergy.setProperty("value", 1.5)
		self.Metal1FermiEnergy.setObjectName( fromUtf8("Metal1FermiEnergy"))
		self.verticalLayout_9.addWidget(self.Metal1FermiEnergy)
		self.label_7 = QtGui.QLabel(self.Metal1Group)
		self.label_7.setObjectName( fromUtf8("label 7"))
		self.verticalLayout_9.addWidget(self.label_7)
		self.Metal1Lattice = QtGui.QDoubleSpinBox(self.Metal1Group)
		self.Metal1Lattice.setDecimals(5)
		self.Metal1Lattice.setObjectName( fromUtf8("Metal1Lattice"))
		self.verticalLayout_9.addWidget(self.Metal1Lattice)
		self.horizontalLayout_3.addLayout(self.verticalLayout_9)
		self.horizontalLayout_4.addWidget(self.Metal1Group)
		self.FerroGroup = QtGui.QGroupBox(self.layoutWidget1)
		self.FerroGroup.setObjectName( fromUtf8("FerroGroup"))
		self.horizontalLayout_5 = QtGui.QHBoxLayout(self.FerroGroup)
		self.horizontalLayout_5.setObjectName( fromUtf8("horizontalLayout_5"))
		self.verticalLayout_11 = QtGui.QVBoxLayout()
		self.verticalLayout_11.setObjectName( fromUtf8("verticalLayout_11"))
		self.label_12 = QtGui.QLabel(self.FerroGroup)
		self.label_12.setObjectName( fromUtf8("label 12"))
		self.verticalLayout_11.addWidget(self.label_12)
		self.FerroDielectricConst = QtGui.QDoubleSpinBox(self.FerroGroup)
		self.FerroDielectricConst.setMaximum (500000.0)
		self.FerroDielectricConst.setProperty("value", 2000.0)
		self.FerroDielectricConst.setObjectName( fromUtf8("FerroDielectricConst"))
		self.verticalLayout_11.addWidget(self.FerroDielectricConst)
		self.label_6 = QtGui.QLabel(self.FerroGroup)
		self.label_6.setObjectName( fromUtf8("label 6"))
		self.verticalLayout_11.addWidget(self.label_6)
		self.FerroBandgap = QtGui.QDoubleSpinBox(self.FerroGroup)
		self.FerroBandgap.setProperty("value", 2.3)
		self.FerroBandgap.setObjectName( fromUtf8("FerroBandgap"))
		self.verticalLayout_11.addWidget(self.FerroBandgap)
		self.label_9 = QtGui.QLabel(self.FerroGroup)
		self.label_9.setObjectName( fromUtf8("label 9"))
		self.verticalLayout_11.addWidget(self.label_9)
		self.FerroPolarization = QtGui.QDoubleSpinBox(self.FerroGroup)
		self.FerroPolarization.setProperty("value", 20.0)
		self.FerroPolarization.setObjectName( fromUtf8("FerroPolarization"))
		self.verticalLayout_11.addWidget(self.FerroPolarization)
		self.label_13 = QtGui.QLabel(self.FerroGroup)
		self.label_13.setObjectName( fromUtf8("label 13"))
		self.verticalLayout_11.addWidget(self.label_13)
		self.FerroThick = QtGui.QDoubleSpinBox(self.FerroGroup)
		self.FerroThick.setProperty("value", 2.0)
		self.FerroThick.setObjectName( fromUtf8("FerroThick"))
		self.verticalLayout_11.addWidget(self.FerroThick)
		self.horizontalLayout_5.addLayout(self.verticalLayout_11)
		self.horizontalLayout_4.addWidget(self.FerroGroup)
		self.Metal2Group = QtGui.QGroupBox(self.layoutWidget1)
		self.Metal2Group.setObjectName( fromUtf8("Metal2Group"))
		self.horizontalLayout_6 = QtGui.QHBoxLayout(self.Metal2Group)
		self.horizontalLayout_6.setObjectName( fromUtf8("horizontalLayout_6"))
		self.verticalLayout_10 = QtGui.QVBoxLayout()
		self.verticalLayout_10.setObjectName( fromUtf8("verticalLayout_10"))
		self.label_11 = QtGui.QLabel(self.Metal2Group)
		self.label_11.setObjectName( fromUtf8("label 11"))
		self.verticalLayout_10.addWidget(self.label_11)
		self.Metal2ScreenLength = QtGui.QDoubleSpinBox(self.Metal2Group)
		self.Metal2ScreenLength.setDecimals(5)
		self.Metal2ScreenLength.setProperty("value", 0.07)
		self.Metal2ScreenLength.setObjectName( fromUtf8("Metal2ScreenLength"))
		self.verticalLayout_10.addWidget(self.Metal2ScreenLength)
		self.label_8 = QtGui.QLabel(self.Metal2Group)
		self.label_8.setObjectName( fromUtf8("label 8"))
		self.verticalLayout_10.addWidget(self.label_8)
		self.Metal2FermiEnergy = QtGui.QDoubleSpinBox(self.Metal2Group)
		self.Metal2FermiEnergy.setProperty("value", 3.5)
		self.Metal2FermiEnergy.setObjectName( fromUtf8("Metal2FermiEnergy"))
		self.verticalLayout_10.addWidget(self.Metal2FermiEnergy)
		self.label_10 = QtGui.QLabel(self.Metal2Group)
		self.label_10.setObjectName( fromUtf8("label 10"))
		self.verticalLayout_10.addWidget(self.label_10)
		self.Metal2Lattice = QtGui.QDoubleSpinBox(self.Metal2Group)
		self.Metal2Lattice.setDecimals(5)
		self.Metal2Lattice.setObjectName( fromUtf8("Metal2Lattice"))
		self.verticalLayout_10.addWidget(self.Metal2Lattice)
		self.horizontalLayout_6.addLayout(self.verticalLayout_10)
		self.horizontalLayout_4.addWidget(self.Metal2Group)
		self.verticalLayout_12.addLayout(self.horizontalLayout_4)
		self.splitter_3 = QtGui.QSplitter(self.splitter_4)
		self.splitter_3.setOrientation(QtCore.Qt.Horizontal)
		self.splitter_3.setObjectName( fromUtf8("splitter 3"))
		self.RunModeling = QtGui.QPushButton(self.splitter_3)
		self.RunModeling.setObjectName( fromUtf8("RunModeling"))
		self.FrontOutput = QtGui.QTextBrowser(self.splitter_3)
		self.FrontOutput.setMaximumSize(QtCore.QSize(16777215, 16777215))
		self.FrontOutput.setObjectName( fromUtf8("FrontOutput"))
		self.verticalLayout_7.addWidget(self.splitter_4)
		self.layoutWidget2 = QtGui.QWidget(self.splitter)
		self.layoutWidget2.setObjectName( fromUtf8("layoutWidget2"))
		self.verticalLayout_14 = QtGui.QVBoxLayout(self.layoutWidget2)
		self.verticalLayout_14.setObjectName( fromUtf8("verticalLayout_14"))
		self.splitter_2 = QtGui.QSplitter(self.layoutWidget2)
		self.splitter_2.setOrientation(QtCore.Qt.Vertical)
		self.splitter_2.setObjectName( fromUtf8("splitter 2"))
		self.verticalWidget_6 = QtGui.QWidget(self.splitter_2)
		self.verticalWidget_6.setObjectName( fromUtf8("verticalWidget_6"))
		self.verticalLayout_8 = QtGui.QVBoxLayout(self.verticalWidget_6)
		self.verticalLayout_8.setObjectName( fromUtf8("verticalLayout_8"))
		self.PotentialMPlay = QtGui.QVBoxLayout()
		self.PotentialMPlay.setObjectName( fromUtf8("PotentialMPlay"))
		self.verticalLayout_8.addLayout(self.PotentialMPlay)
		self.PotentialOutput = QtGui.QTextBrowser(self.splitter_2)
		self.PotentialOutput.setMaximumSize(QtCore.QSize(16777215, 16777215))
		self.PotentialOutput.setObjectName( fromUtf8("PotentialOutput"))
		self.verticalLayout_14.addWidget(self.splitter_2)
		self.widget = QtGui.QWidget(self.splitter_5)
		self.widget.setObjectName( fromUtf8("widget"))
		self.verticalLayout_2 = QtGui.QVBoxLayout(self.widget)
		self.verticalLayout_2.setObjectName( fromUtf8("verticalLayout_2"))
		self.ModelingProgress = QtGui.QProgressBar(self.widget)
		self.ModelingProgress.setProperty("value", 24)
		self.ModelingProgress.setObjectName( fromUtf8("ModelingProgress"))
		self.verticalLayout_2.addWidget(self.ModelingProgress)
		self.horizontalLayout_7 = QtGui.QHBoxLayout()
		self.horizontalLayout_7.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
		self.horizontalLayout_7.setObjectName( fromUtf8("horizontalLayout_7"))
		self.label_16 = QtGui.QLabel(self.widget)
		self.label_16.setObjectName( fromUtf8("label 16"))
		self.horizontalLayout_7.addWidget(self.label_16)
		self.label_15 = QtGui.QLabel(self.widget)
		self.label_15.setObjectName( fromUtf8("label 15"))
		self.horizontalLayout_7.addWidget(self.label_15)
		self.verticalLayout_2.addLayout(self.horizontalLayout_7)
		self.gridLayout.addWidget(self.splitter_5 , 0, 0, 1, 1)
		MainWindow.setCentralWidget(self.centralwidget)
		self.statusbar = QtGui.QStatusBar(MainWindow)
		self.statusbar.setObjectName( fromUtf8("statusbar"))
		MainWindow.setStatusBar(self.statusbar)
		self.retranslateUi(MainWindow)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)

	def retranslateUi(self, MainWindow):
		MainWindow.setWindowTitle( translate("MainWindow", "MainWindow", None))
		self.groupBox.setTitle( translate("MainWindow", "Materials", None))
		self.label.setText( translate("MainWindow", "Metal 1", None))
		self.label_2.setText( translate("MainWindow", "Ferroelectric", None))
		self.label_3.setText( translate("MainWindow", "Metal 2", None))
		self.label_14.setText( translate("MainWindow", "If you know the metalscreening lengths , enter them below and leave lattice constantas 0.\n Otherwise , leave them as 0 and enter the lattice constant , and the program will calculate an approx. screening length.", None))
		self.Metal1Group.setTitle( translate("MainWindow", "Metal 1", None))
		self.label_4.setText( translate("MainWindow", "Screening Length (nm)", None))
		self.label_5.setText( translate("MainWindow", "Fermi Energy (eV)",None))
		self.label_7.setText( translate("MainWindow", "Lattice Constant (Angstrom)", None))
		self.FerroGroup.setTitle( translate("MainWindow", "Ferroelectric",None))
		self.label_12.setText( translate("MainWindow", "Dielectric Constant (E f / E 0 )", None))
		self.label_6.setText( translate("MainWindow", "Bandgap (eV)", None))
		self.label_9.setText( translate("MainWindow", "Polarization (micro C/cm^2)", None))
		self.label_13.setText( translate("MainWindow", "Thickness (nm)", None))
		self.Metal2Group.setTitle( translate("MainWindow", "Metal 2", None))
		self.label_11.setText( translate("MainWindow", "Screening Length (nm)", None))
		self.label_8.setText( translate("MainWindow", "Fermi Energy (eV)",None))
		self.label_10.setText( translate("MainWindow", "Lattice Constant (Angstrom)", None))
		self.RunModeling.setText( translate("MainWindow", "Model This!\n"
		"\n"
		"Output values\n"
		"will appear to the\n"
		" right--->\n"
		"\n"
		"Energy band is plotted\n"
		"at far right.\n"
		"\n"
		"Further plots will be found\n"
		" on other tabs.", None))
		self.label_16.setText( translate("MainWindow", "May 2016, sap1951@rit.edu, (585) 236-9510", None))
		self.label_15.setText( translate("MainWindow", "Created by SpencerPringle for Rochester Institute of Technology , 1 Lomb Drive ,Rochester , NY 14623", None))
