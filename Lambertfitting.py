 # -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 15:51:07 2014

@author: aspitarl
"""

#from __future__ import division
#from pylab import *
#import Tkinter
from Tkinter import *
#from Tkinter import Tk, Frame, IntVar
#from ttk import  Style, Entry
#import sys
import tkFileDialog
import matplotlib
#matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#from matplotlib.figure import F,igure
#import matplotlib.pyplot as pltlib

from scipy.special import lambertw
from scipy import zeros
import numpy as np
#import math
#from lmfit import minimize, Parameters,Parameter, report_fit, Minimizer
from lmfit import *


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def lamb(x, i0,nvt,R,Voff): 
    return i0*((nvt/(i0*R))*lambertw((i0*R/nvt)*np.exp(((x-Voff)+i0*R)/nvt)).real - 1)

def fnc2min(params,x,data):
    i0 = params['i0'].value
    nvt = params['nvt'].value
    R = params['R'].value
    Voff = params['Voff'].value
    return lamb(x,i0,nvt,R,Voff)-data
    

  
    
    
class Example(Frame):
    
    def __init__(self,parent):
        Frame.__init__(self,parent,background="white")
        self.parent = parent
        self.initUI()
    
    def initUI(self):
        self.parent.title("Curve Fitting")
        
        self.file_path = ""
        #declare global variables
#        self.popt = zeros(3)
        self.popt = Parameters()
        self.results = minimizer.MinimizerResult()
        self.d = zeros((2,2))
        
        self.vdata = zeros(2)
        self.test = zeros((2,2))
        self.fit = zeros((2,2))        
        self.usei0=IntVar()
        self.usenvt=IntVar()
        self.useR=IntVar()
        self.useVoff=IntVar()
        self.usetest=IntVar()
        self.newfit=IntVar()
        self.log=IntVar()
        self.lin=IntVar()
        
        cblabel= Label(self, text="fit Value?")
        cblabel.grid(row=0, column=0)


                
        
        #setup plot space
    #    self.canvasFig=pltlib.figure(1)
        self.Fig = matplotlib.figure.Figure(figsize=(5,4),dpi=100)
        self.FigSubPlot = self.Fig.add_subplot(111)
        #add plot
        self.line1, = self.FigSubPlot.plot(self.d[:,0],self.d[:,1],'r-')
        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.Fig, master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.canvas.get_tk_widget().grid(column=7,row=0, rowspan =5,columnspan=2)
        self.update()#??
              
              
        labelfitin = Label(self, text="fit input", height=4)
        labelfitin.grid(row=1, column=0)
#        i0

        self.entryi0 = Entry(self)
        self.entryi0.grid(row=1, column=1, sticky=W+E)
        self.entryi0.insert(0,1e-12)                
        self.cbi0= Checkbutton(self,text="i0",variable=self.usei0)
        self.cbi0.select()
        self.cbi0.grid(row=0,column=1)
        
#        nvt
        self.entrynvt = Entry(self)
        self.entrynvt.grid(row=1, column=2, sticky=W+E)
        self.entrynvt.insert(0,.025)        
        self.cbnvt= Checkbutton(self,text="nVt",variable=self.usenvt)
        self.cbnvt.select()
        self.cbnvt.grid(row=0,column=2)


#       R
        self.entryR = Entry(self)
        self.entryR.grid(row=1, column=3, sticky=W+E)
        self.entryR.insert(0,"%g" %(1e+9))
        self.cbR= Checkbutton(self,text="R",variable=self.useR)
        self.cbR.select()
        self.cbR.grid(row=0,column=3)
#       Voff
        self.entryVoff = Entry(self)
        self.entryVoff.grid(row=1, column=4, sticky=W+E)
        self.entryVoff.insert(0,"%g" %(0))
        self.cbVoff= Checkbutton(self,text="Voff",variable=self.useVoff)
        self.cbVoff.select()
        self.cbVoff.grid(row=0,column=4)
#        test function before fitting 
        self.cbtest= Checkbutton(self,text="test function",variable=self.usetest)
        self.cbtest.select()
        self.cbtest.grid(row=4,column=1)
        self.testline, = self.FigSubPlot.plot(self.test[:,0],self.test[:,1],'r-')
        
        self.fitline, = self.FigSubPlot.plot(self.fit[:,0],self.fit[:,1],'r-')
        
    
        labelfitout = Label(self, text="Fit output", height=4)
        labelfitout.grid(row=2, column=0)
        
        self.fiti0 = Text(self, width=12, height =1)
        self.fiti0.grid(row=2,column=1)
        self.fitnvt = Text(self, width=12, height =1)
        self.fitnvt.grid(row=2,column=2)
        self.fitR = Text(self, width=12, height =1)
        self.fitR.grid(row=2,column=3)
        self.fitVoff = Text(self, width=12, height =1)
        self.fitVoff.grid(row=2,column=4)
#       Input for bounds for curve fitting

        labelminV=Label(self, text="minV", height=4)
        labelminV.grid(row=3, column=0)
        labelmaxV=Label(self, text="maxV", height=4)
        labelmaxV.grid(row=3, column=2)       
        self.entryminV = Entry(self)
        self.entryminV.grid(row=3, column=1, sticky=W+E)
        self.entryminV.insert(0,0.25)
        self.entrymaxV = Entry(self)
        self.entrymaxV.grid(row=3, column=3, sticky=W+E)     
        self.entrymaxV.insert(0,0.5)
       
        quitbutton= Button(self,text="quit", command = self.parent.destroy)
        quitbutton.grid(row=4,column=0)
        
        refreshbutton= Button(self,text="refresh", command = self.refreshFigure)
        refreshbutton.grid(row=4,column=2)
        
        loadbutton= Button(self,text="loadfile", command = self.loadfile)
        loadbutton.grid(row=4,column=3)

        fitbutton= Button(self,text="Fit!", command = self.dofit)
        fitbutton.grid(row=5,column=2)
        
        transbutton= Button(self,text="Trans parameters", command = self.transfer)
        transbutton.grid(row=5,column=3)
        
        
        labelcalcnvt=Label(self, text="T,n:", height=4)
        labelcalcnvt.grid(row=6, column=0)
        self.entryT = Entry(self)
        self.entryT.grid(row=6, column=1, sticky=W+E)
        self.entryT.insert(0,300)
        self.entryn = Entry(self)
        self.entryn.grid(row=6, column=2, sticky=W+E)     
        self.entryn.insert(0,1)        
        calcnvtbutton= Button(self,text="calc nvt", command = self.calcnvt)
        calcnvtbutton.grid(row=6,column=3)
        
        self.loggraph= Checkbutton(self,text="log",variable=self.log)
        self.loggraph.deselect()
        self.loggraph.grid(row=6,column=6)
        
        self.lingraph= Checkbutton(self,text="linear",variable=self.lin)
        self.lingraph.select()
        self.lingraph.grid(row=6,column=7)
        
        addavgbutton= Button(self,text="add average", command = self.addaverage)
        addavgbutton.grid(row=6,column=8)
        
        printbutton= Button(self,text="print array", command = self.printarray)
        printbutton.grid(row=6,column=9)
        
        cutbutton= Button(self,text="cut array", command = self.cutarray)
        cutbutton.grid(row=6,column=10)
        
        cutbutton= Button(self,text="save", command = self.save)
        cutbutton.grid(row=7,column=10)        
        
        self.pack()
        
    def refreshFigure(self):
        if len(self.FigSubPlot.lines) > 1 : #if there is a plot remove it
            for i in range(len(self.FigSubPlot.lines)):
                self.FigSubPlot.lines[-1].remove()
        self.line1, = self.FigSubPlot.plot(self.d[:,0],self.d[:,1],'o',color='g') #add new plot
#        test function is no longer nesssecary as can fit with zero params
        testvdata = np.arange(self.d[:,0].min(),self.d[:,0].max(),(self.d[1,0]-self.d[0,0]))
        if (self.usetest.get() == 1):                    
            self.linetest, = self.FigSubPlot.plot(testvdata,lamb(testvdata,float(self.entryi0.get()),float(self.entrynvt.get()),float(self.entryR.get()),float(self.entryVoff.get())),'r-',color='b')
            
        if (self.newfit == 1):
            self.fitline, = self.FigSubPlot.plot(testvdata,lamb(testvdata,self.results.params['i0'].value,self.results.params['nvt'].value,self.results.params['R'].value,self.results.params['Voff'].value),color = 'r')
            self.fiti0.delete("1.0",END)
            self.fiti0.insert(INSERT, "%g" %(self.results.params['i0'].value))
            self.fitnvt.delete("1.0",END)
            self.fitnvt.insert(INSERT, "%g" %(self.results.params['nvt'].value))
            self.fitR.delete("1.0",END)
            self.fitR.insert(INSERT, "%g" %(self.results.params['R'].value))
            self.fitVoff.delete("1.0",END)
            self.fitVoff.insert(INSERT, "%g" %(self.results.params['Voff'].value))
        ax = self.FigSubPlot.axes
        if(self.log.get() ==1):
            ax.set_xlim(self.d[:,0].min(), self.d[:,0].max())   #set up axes
            ax.set_ylim(abs(self.d[:,1]).min(), self.d[:,1].max())
            self.FigSubPlot.set_yscale('log')
            self.lingraph.deselect()
        if(self.lin.get() == 1):
            ax.set_xlim(self.d[:,0].min(), self.d[:,0].max())   #set up axes
            ax.set_ylim(self.d[:,1].min(), self.d[:,1].max())
            self.FigSubPlot.set_yscale('linear')
        
        self.FigSubPlot.axvline(self.entryminV.get())  #add curve fit bound lines
        self.FigSubPlot.axvline(self.entrymaxV.get())
        
        self.canvas.draw()
        
        #http://tkinter.unpythonic.net/wiki/tkFileDialog
        #https://mail.python.org/pipermail/tutor/2012-February/088267.html
    def loadfile(self):
        self.file_path = tkFileDialog.askopenfilename(initialdir='C:/Users/')
        self.d = np.loadtxt(self.file_path, delimiter="\t",skiprows=2, usecols= (0,1))
        self.refreshFigure()
        
    def dofit(self):
        minVidx = find_nearest(self.d[:,0],float(self.entryminV.get()))
        maxVidx = find_nearest(self.d[:,0],float(self.entrymaxV.get()))
        self.vdata = self.d[minVidx:maxVidx,0]
        idata = self.d[minVidx:maxVidx,1]
        

            
        self.popt['i0'] = Parameter(value=float(self.entryi0.get()),vary=self.usei0.get())
        self.popt['nvt'] = Parameter(value=float(self.entrynvt.get()), vary =self.usenvt.get() )
        self.popt['R'] = Parameter(value=float(self.entryR.get()), vary = self.useR.get())
        self.popt['Voff'] = Parameter(value=float(self.entryVoff.get()), vary = self.useVoff.get())
                
        self.results= minimize(fnc2min, self.popt, args=(self.vdata,idata))
        print self.results
        self.newfit=1
        self.refreshFigure()
    def transfer(self):
        if (self.usei0.get()==0):
            self.entryi0.delete(0,END)
            self.entryi0.insert(0, "%g" %(self.results.params['i0'].value))
        if (self.usenvt.get()==0):
            self.entrynvt.delete(0,END)
            self.entrynvt.insert(0, "%g" %(self.results.params['nvt'].value))
        if (self.useR.get()==0):
            self.entryR.delete(0,END)
            self.entryR.insert(0, "%g" %(self.results.params['R'].value))
        if (self.useVoff.get()==0):
            self.entryVoff.delete(0,END)
            self.entryVoff.insert(0, "%g" %(self.results.params['Voff'].value))
    def calcnvt(self):
        self.entrynvt.delete(0,END)
        self.entrynvt.insert(0, (.026)*(float(self.entryT.get())/300)*float(self.entryn.get()))

    def addaverage(self):
        minVidx = find_nearest(self.d[:,0],float(self.entryminV.get()))
        maxVidx = find_nearest(self.d[:,0],float(self.entrymaxV.get()))
        drange = self.d[minVidx:maxVidx,1]
        dmean = np.mean(drange)
        print dmean
        self.d= self.d - dmean
    def printarray(self):
        print "printing array"
        print self.d
    def save(self):
        filename = tkFileDialog.asksaveasfilename(initialdir='C:/Users/')        
        vdata = np.arange(self.d[:,0].min(),self.d[:,0].max(),(self.d[1,0]-self.d[0,0]))
        Idata = lamb(vdata,self.results.params['i0'].value,self.results.params['nvt'].value,self.results.params['R'].value,self.results.params['Voff'].value)
        data = np.dstack((Idata, vdata))#round about way since i don't understand
        data = np.squeeze(data)
        
        Tval = float(self.entryT.get())
        nval = self.results.params['nvt'].value/(0.026*(Tval/300))
        
        headerstr = "Filepath: " + self.file_path
        headerstr += "\nI0: " + str(self.results.params['i0'].value)
        
        headerstr += " nVT: " + str(self.results.params['nvt'].value)
        headerstr += " T: " + self.entryT.get()
        headerstr += " n: " + str(nval)
        
        headerstr += " R: " + str(self.results.params['R'].value)        
        headerstr += " Voff: " + str(self.results.params['Voff'].value)  
        
        np.savetxt(filename, data, header = headerstr)
        
    def cutarray(self):
        exit =0
        i=0
        while((exit ==0) and (i<len(self.d[:,0]))):
            if(self.d[i,0]>self.d[i+1,0]):
                self.d= np.delete(self.d,np.arange(i,len(self.d[:,0])),0)
#                del self.d[i:,1]
                print self.d
                exit==1
            i=i+1
            
root = Tk()
#root.geometry("250x250+300+300")

w=Example(root)

root.mainloop()


#sys.exit(0)
