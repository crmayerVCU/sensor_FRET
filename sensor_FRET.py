"""
Author: Carl Mayer, crmayer@vcu.edu

See ???publication??? for details regarding implentation

For this script you need to acquire images of your FRET sensor at two distinct 
excitation wavelengths, eg. for a Cerulean-Venus FRET pair you would need an 
image of the same cell/region taken with both 405nm and 445nm lasers to calculate 
the FRET efficiency.

IMPORTANT:
The input to this script should be a 4 channel image structured in the following order:
[Donor @ ex1, Acceptor @ ex1, Donor @ ex2, Acceptor @ ex2]

The the donor and acceptor channels can be determined using linear unmixing, or if
a spectral detector is unavailable (eg. widefield or airyscan data where physical 
filters are used) the closest available bandpass channels can be used as long as the
single fluorophore basis vectors are defined (see publication above for details, in 
the case of already unmixed data the basis vectors will be [1,0],[0,1] for both 
excitations).
 
If needed, additional non-FRET related fluorophores (for supplemental labeling) can be 
extra channels after the first 4 without affecting the code.

Any ROIs called out in the ROI Manager window will also be analyzed and saved in a table. 
ROI based analysis is preferable due to the nature of the Cauchy distributed noise in the 
calculated FRET values, as averaging fluorophore intensity and then calculating FRET is 
more accurate than simply averaging the FRET determined for each individual pixel.

Results are saved in a new directory labeled “/FileBaseName_Results/“ and contain calculated 
images, a table with ROI results, and a text file stating the calibration values used in 
the calculations.
"""

from ij import IJ, ImagePlus, ImageStack
from ij.process import FloatProcessor, LUT
from ij.process import ImageStatistics as IS
from ij.plugin import ImageCalculator, Duplicator
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.io import FileSaver
from ij.gui import GenericDialog  
import sys
import os
import math
import re
import csv
 

"""
Read In Default Data
"""

#Calibration Data: adjust for fluorophores and excitation wavelengths used

#Fluorophore Names
donor='Cerulean' 
acceptor='Venus'
#Quantum Yields for each fluorophore
Qd=0.62  
Qa=0.57 
#excitation laser wavelength names
ex1='405'  
ex2='445' 
#Ratio of extinction coefficients at each excitation wavelength, eD/eA
ex1_ratio= 17.105 
ex2_ratio= 2.654  
gamma=ex2_ratio/ex1_ratio #Calibration parameter
#Donor and acceptor basis vectors w/ respect to measured bandpass channels at each excitation wavelength
#formatted Dx,Dy,Ax,Ay
B1string='4882.271,1266.419,80.552,71.149'		
B2string='3674.368,1305.003,382.355,374.865'
#Intensity offset to subtract (Airyscan images have a baseline shift parameter that adds 10000 intensity units to every pixel)
offset=10000


#setup GUI for Calibration Data Input
def calibOptions():  
  gd = GenericDialog("User Input")  
  gd.addStringField("Donor Fluorophore", donor) 
  gd.addStringField("Acceptor Fluorophore", acceptor)  
  gd.addNumericField("Donor Quantum Yield", Qd, 3)
  gd.addNumericField("Acceptor Quantum Yield", Qa, 3) 
  gd.addStringField("Excitation Wavelength 1", ex1) 
  gd.addStringField("Excitation Wavelength 2", ex2)
  gd.addNumericField("Extinction Coef. Ratio 1", ex1_ratio, 3)
  gd.addNumericField("Extinction Coef. Ratio 2", ex2_ratio, 3)
  gd.addStringField("Basis Vectors 1", B1string) 
  gd.addStringField("Basis Vectors 2", B2string) 
  gd.addNumericField("Gamma", gamma, 3) 
  gd.addNumericField("Offset", offset, 3) 
  gd.showDialog()  
  #  
  if gd.wasCanceled():    
    sys.exit('Cancelled')
    return  
    
  # Read out the options  
  d = gd.getNextString()  
  a = gd.getNextString()  
  qd = gd.getNextNumber() 
  qa = gd.getNextNumber()  
  e1 = gd.getNextString()  
  e2 = gd.getNextString() 
  er1 = gd.getNextNumber() 
  er2 = gd.getNextNumber()
  b1 = gd.getNextString()
  b2 = gd.getNextString() 
  g = gd.getNextNumber()
  o = gd.getNextNumber() 
  
  return d,a,qd,qa,e1,e2,er1,er2,b1,b2,g,o
  
options = calibOptions()  
if options is not None:  
  donor,acceptor,Qd,Qa,ex1,ex2,ex1_ratio,ex2_ratio,B1string,B2string,gamma,offset = options

B1=[float(i) for i in options[8].split(',')]
B2=[float(i) for i in options[9].split(',')]
B1coef=[(B1[0]**2+B1[1]**2)**0.5*B1[3]/(B1[0]*B1[3]-B1[2]*B1[1]), 
		(B1[0]**2+B1[1]**2)**0.5*B1[2]/(B1[0]*B1[3]-B1[2]*B1[1]),
		(B1[2]**2+B1[3]**2)**0.5*B1[0]/(B1[0]*B1[3]-B1[2]*B1[1]), 
		(B1[2]**2+B1[3]**2)**0.5*B1[1]/(B1[0]*B1[3]-B1[2]*B1[1])]
B2coef=[(B2[0]**2+B2[1]**2)**0.5*B2[3]/(B2[0]*B2[3]-B2[2]*B2[1]), 
		(B2[0]**2+B2[1]**2)**0.5*B2[2]/(B2[0]*B2[3]-B2[2]*B2[1]),
		(B2[2]**2+B2[3]**2)**0.5*B2[0]/(B2[0]*B2[3]-B2[2]*B2[1]), 
		(B2[2]**2+B2[3]**2)**0.5*B2[1]/(B2[0]*B2[3]-B2[2]*B2[1])]
print(B1coef)
print(B2coef)

#read in already opened image and read img dimensions
imp=IJ.getImage()
folder=imp.getOriginalFileInfo().directory
fn=imp.getOriginalFileInfo().fileName
h=imp.height
w=imp.width
print(folder)

#pull out base file name
fnb=fn.rpartition('.')[0]
print(fnb)
#make a directory to store results in
resPath=os.path.join(folder,fnb+'_Results')
try:
	os.mkdir(resPath)
except:
	print('Results already exist, overwriting...')
print(resPath)

#dumps a record of the calibration data used 
#into a text file in results folder
f = open(os.path.join(resPath,fnb+'_settings.csv'), 'wb')
writer = csv.writer(f)

writer.writerow(["Donor Fluorophore", donor]) 
writer.writerow(["Acceptor Fluorophore", acceptor])  
writer.writerow(["Donor Quantum Yield", Qd])
writer.writerow(["Acceptor Quantum Yield", Qa]) 
writer.writerow(["Excitation Wavelength 1", ex1]) 
writer.writerow(["Excitation Wavelength 2", ex2])
writer.writerow(["Extinction Coef. Ratio 1", ex1_ratio])
writer.writerow(["Extinction Coef. Ratio 2", ex2_ratio])
writer.writerow(["Basis Vector 1", B1string])
writer.writerow(["Basis Vector 2", B2string]) 
writer.writerow(["Gamma", gamma]) 
writer.writerow(["Offset", offset])
f.close()


"""
Internally Used Functions
"""
#gets statistical values (mean median etc.) for an image
def imstats(imp):
	ip=imp.getProcessor()
	options=IS.MEAN | IS.MEDIAN | IS.MIN_MAX
	stats=IS.getStatistics(ip,options,imp.getCalibration())
	return stats

#calculates various FRET parameters for a single observation
def FRETcalc(X1,X2,Y1,Y2):
	#Remove offset and convert bandpass channel intensities to fluorophore intensities using fluorophore basis vectors
	X1=X1-offset
	X2=X2-offset
	Y1=Y1-offset
	Y2=Y2-offset
	D1=B1coef[0]*X1-B1coef[1]*Y1
	D2=B2coef[0]*X2-B2coef[1]*Y2
	A1=B1coef[2]*Y1-B1coef[3]*X1
	A2=B2coef[2]*Y2-B2coef[3]*X2
	
	#Average Intensity
	I=(D1+D2+A1+A2)/4

	#Qualitative FRET Index Values for each excitation
	Ind1=A1/D1
	Ind2=A2/D2

	#Intermediate values to get to the FRET Efficiency
	Alpha=D2/D1
	Beta=A2-(Alpha*A1)
	Adir1=Beta/(Alpha*(gamma**-1-1))
	Adir2=Beta/(1-gamma)
	eDeA1=(D1/Qd+(A1-Adir1)/Qa)/(Adir1/Qa)
	eDeA2=(D2/Qd+(A2-Adir2)/Qa)/(Adir2/Qa)

	#negative extinction coefficients are ignored
	
	
	if eDeA2<0 or eDeA1<0:
		S=float('nan')
		E_exD=float('nan')
		E_exA=float('nan')
	else:
		#Stoichemetry ([D] to [A] ratio on a log scale)
		S=math.log(eDeA2/ex2_ratio,10)
		#FRET efficiency assuming excess Donor
		E_exD=((A2-Adir2)/Adir2)/ex2_ratio
		#FRET Efficiency assuming excess Acceptor
		E_exA=(((D2*Qa)/((A2-Adir2)*Qd))+1)**-1

	

	#Best estimate of the FRET Efficiency
	E=max(E_exD,E_exA)

	return(D1,D2,A1,A2,I,Ind1,Ind2,Alpha,Beta,Adir1,Adir2,
			eDeA1,eDeA2,S,E_exD,E_exA,E)

#converts the ImagePlus object to a pixel array for analysis
def IMPtoPIX(imp):
	imp_ip = imp.getProcessor().convertToFloat() # as a copy
	imp_pix = imp_ip.getPixels()
	return imp_pix

#converts the pixel array back to an ImagePlus object for display and native imageJ functions
def PIXtoIMP (imp_pix,title,height,width):
	ip = FloatProcessor(width, height, imp_pix, None)
	imp = ImagePlus(title, ip)
	return imp

#saves images to results folder
def impSave(imp):
	imp.setProperty('DonorFluorophore:','mtfp')
	fs=FileSaver(imp)
	fs.saveAsTiff(os.path.join(resPath, fnb+'_'+imp.title+'.tif'))
	return

#copys columns from one table to another for ROI analysis
def copyColumn(RTfrom,RTto,fromColName,toColName):
	for i in range(RTfrom.size()):
		RTto.setValue(toColName,i,RTfrom.getValue(fromColName,i))
	return(RTto)

#performs the FRET calculations for each ROI and populates the results table
def roiFRET(RT):
	for i in range(RT.size()):
		X1roi=RT.getValue(X1.title,i)
		X2roi=RT.getValue(X2.title,i)
		Y1roi=RT.getValue(Y1.title,i)
		Y2roi=RT.getValue(Y2.title,i)
		vals=FRETcalc(X1roi,X2roi,Y1roi,Y2roi)
		RT.setValue(D1.title,i,vals[0])
		RT.setValue(D2.title,i,vals[1])
		RT.setValue(A1.title,i,vals[2])
		RT.setValue(A2.title,i,vals[3])
		RT.setValue(I.title,i,vals[4])
		RT.setValue(Ind1.title,i,vals[5])
		RT.setValue(Ind2.title,i,vals[6])
		RT.setValue(Alpha.title,i,vals[7])
		RT.setValue(Beta.title,i,vals[8])
		RT.setValue(Adir1.title,i,vals[9])
		RT.setValue(Adir2.title,i,vals[10])
		RT.setValue(eDeA1.title,i,vals[11])
		RT.setValue(eDeA2.title,i,vals[12])
		RT.setValue(S.title,i,vals[13])
		RT.setValue(E_exD.title,i,vals[14])
		RT.setValue(E_exA.title,i,vals[15])
		RT.setValue(E.title,i,vals[16])
	return


"""
Analysis
"""
#Split out each FRET channel(X and Y) for both excitation frequencies
#Images must be formatted [X1,Y1,X2,Y2,N1,N2,N3...] where any N channel
#can be used to quantify autofluorescence or Fluorescent labels unrelated to the FRET
X1=Duplicator().run(imp,1,1)
X1.title='Ch-X '+ex1+'nm'
X1_pix=IMPtoPIX(X1)

X2=Duplicator().run(imp,3,3)
X2.title='Ch-X '+ex2+'nm'
X2_pix=IMPtoPIX(X2)

Y1=Duplicator().run(imp,2,2)
Y1.title='Ch-Y '+ex1+'nm'
Y1_pix=IMPtoPIX(Y1)

Y2=Duplicator().run(imp,4,4)
Y2.title='Ch-Y '+ex2+'nm'
Y2_pix=IMPtoPIX(Y2)
		
#Calculate FRET related values from each pixel and imput them into
#empty pixel lists to fill for calculated images
D1_pix=[]
D2_pix=[]
A1_pix=[]
A2_pix=[]
I_pix=[]
Ind1_pix=[]
Ind2_pix=[]
Alpha_pix=[]
Beta_pix=[]
Adir1_pix=[]
Adir2_pix=[]
eDeA1_pix=[]
eDeA2_pix=[]
S_pix=[]
E_exD_pix=[]
E_exA_pix=[]
E_pix=[]
NaN=[float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')]
vals=NaN
count=0
for i in range(len(X1_pix)):
	IJ.showProgress(i, len(X1_pix)+1)
	vals=FRETcalc(X1_pix[i],X2_pix[i],Y1_pix[i],Y2_pix[i])
	if vals[-1]==float('nan'):
		count=count+1
			
		'''
		try:
			vals=FRETcalc(X1_pix[i],X2_pix[i],Y1_pix[i],Y2_pix[i])
		except:
			vals=NaN
			count=count+1
		'''
			
	D1_pix.append(vals[0])
	D2_pix.append(vals[1])
	A1_pix.append(vals[2])
	A2_pix.append(vals[3])
	I_pix.append(vals[4])
	Ind1_pix.append(vals[5])
	Ind2_pix.append(vals[6])
	Alpha_pix.append(vals[7])
	Beta_pix.append(vals[8])
	Adir1_pix.append(vals[9])
	Adir2_pix.append(vals[10])
	eDeA1_pix.append(vals[11])
	eDeA2_pix.append(vals[12])
	S_pix.append(vals[13])
	E_exD_pix.append(vals[14])
	E_exA_pix.append(vals[15])
	E_pix.append(vals[16])

print('NaN Pixels: '+str(count))
#convert each generated pixel list back to an ImagePlus object
D1=PIXtoIMP(D1_pix,donor+' '+ex1+'nm',h,w)
D2=PIXtoIMP(D2_pix,donor+' '+ex2+'nm',h,w)
A1=PIXtoIMP(A1_pix,acceptor+' '+ex1+'nm',h,w)
A2=PIXtoIMP(A2_pix,acceptor+' '+ex2+'nm',h,w)
I=PIXtoIMP(I_pix,'Intensity',h,w)
Ind1=PIXtoIMP(Ind1_pix,'FRET Index '+ex1+'nm',h,w)
Ind2=PIXtoIMP(Ind2_pix,'FRET Index '+ex2+'nm',h,w)
Alpha=PIXtoIMP(Alpha_pix,'Alpha',h,w)
Beta=PIXtoIMP(Beta_pix,'Beta',h,w)
Adir1=PIXtoIMP(Adir1_pix,'Acceptor Direct Excitation '+ex1+'nm',h,w)
Adir2=PIXtoIMP(Adir2_pix,'Acceptor Direct Excitation '+ex2+'nm',h,w)
eDeA1=PIXtoIMP(eDeA1_pix,'Extinction Coef. Ratio '+ex1+'nm',h,w)
eDeA2=PIXtoIMP(eDeA2_pix,'Extinction Coef. Ratio '+ex2+'nm',h,w)
S=PIXtoIMP(S_pix,'Stoicheometry (log([D]_[A])',h,w)
E_exD=PIXtoIMP(E_exD_pix,'FRET Efficiency (excess Donor)',h,w)
E_exA=PIXtoIMP(E_exA_pix,'FRET Efficiency (excess Acceptor)',h,w)
E=PIXtoIMP(E_pix,'FRET Efficiency',h,w)

#sets up gui for you to choose what data to save
def saveOptions():  
  gd = GenericDialog("Save Data")
  gd.addCheckbox(D1.title, False)  
  gd.addCheckbox(D2.title, False)
  gd.addCheckbox(A1.title, False)  
  gd.addCheckbox(A2.title, False)   
  gd.addCheckbox(I.title, False)  
  gd.addCheckbox(Ind1.title, True)  
  gd.addCheckbox(Ind2.title, True)  
  gd.addCheckbox(Alpha.title, False)  
  gd.addCheckbox(Beta.title, False)  
  gd.addCheckbox(Adir1.title, False)  
  gd.addCheckbox(Adir2.title, False)  
  gd.addCheckbox(eDeA1.title, False)  
  gd.addCheckbox(eDeA2.title, False)  
  gd.addCheckbox(S.title, True)  
  gd.addCheckbox(E_exD.title, False)  
  gd.addCheckbox(E_exA.title, False)  
  gd.addCheckbox(E.title, True)  
  gd.showDialog()  
  #  
  if gd.wasCanceled():  
    sys.exit('Cancelled')
    return  
  # Read out the options  
  saveList=[]
  if gd.getNextBoolean():
  	saveList.append(D1.title)
  if gd.getNextBoolean():
  	saveList.append(D2.title)
  if gd.getNextBoolean():
  	saveList.append(A1.title)
  if gd.getNextBoolean():
  	saveList.append(A2.title)
  if gd.getNextBoolean():
  	saveList.append(I.title)
  if gd.getNextBoolean():
  	saveList.append(Ind1.title)
  if gd.getNextBoolean():
  	saveList.append(Ind2.title)
  if gd.getNextBoolean():
  	saveList.append(Alpha.title)
  if gd.getNextBoolean():
  	saveList.append(Beta.title)
  if gd.getNextBoolean():
  	saveList.append(Adir1.title)
  if gd.getNextBoolean():
  	saveList.append(Adir2.title)
  if gd.getNextBoolean():
  	saveList.append(eDeA1.title)
  if gd.getNextBoolean():
  	saveList.append(eDeA2.title)
  if gd.getNextBoolean():
  	saveList.append(S.title)
  if gd.getNextBoolean():
  	saveList.append(E_exD.title)
  if gd.getNextBoolean():
  	saveList.append(E_exA.title)
  if gd.getNextBoolean():
  	saveList.append(E.title)
  	
  return saveList
  
sl = saveOptions()  

"""
#you can manually define what images to save for batch processing,
#just remove triple paranthesis, comment out the sl above and 
#comment out anything in the sl variable below to exclude those 
#images

sl=[
#	D1.title,
#	D2.title,
#	A1.title,
#	A2.title,
#	I.title,
	Ind1.title,
	Ind2.title,
#	Alpha.title,
#	Beta.title,
#	Adir1.title,
#	Adir2.title,
#	eDeA1.title,
#	eDeA2.title,
	S.title,
	E_exD.title,
	E_exA.title,
	E.title
	]
"""

def plotAndSave(imp,rng=None,cmap='phase'):
	
	IJ.run(imp,cmap,"")
	if rng != None:
		imp.setDisplayRange(rng[0],rng[1])
		IJ.run(imp, "Calibration Bar...", "location=[Lower Left] fill=Black label=White number=3 decimal=1 font=12 zoom=1.3 overlay");
	imp.show()
	impSave(imp)
	return

if D1.title in sl:	
	plotAndSave(D1,cmap='Cyan')

if D2.title in sl:	
	plotAndSave(D2,cmap='Cyan')

if A1.title in sl:	
	plotAndSave(A1,cmap='Yellow')
	
if A2.title in sl:	
	plotAndSave(A2,cmap='Yellow')
	
if I.title in sl:	
	plotAndSave(I,cmap='Grays')

if Ind1.title in sl:	
	plotAndSave(Ind1,rng=[0.00,3.00])
	
if Ind2.title in sl:	
	plotAndSave(Ind2,rng=[0.00,3.00])

if Alpha.title in sl:	
	plotAndSave(Alpha,cmap='Grays')

if Beta.title in sl:	
	plotAndSave(Beta,cmap='Grays')

if Adir1.title in sl:	
	plotAndSave(Adir1,cmap='Yellow')

if Adir2.title in sl:	
	plotAndSave(Adir2,cmap='Yellow')

if eDeA1.title in sl:	
	plotAndSave(eDeA1)

if eDeA2.title in sl:	
	plotAndSave(eDeA2)

if S.title in sl:	
	plotAndSave(S,rng=[-1.00,1.00])

if E_exD.title in sl:	
	plotAndSave(E_exD,rng=[0.00,1.00])

if E_exA.title in sl:	
	plotAndSave(E_exA,rng=[0.00,1.00])

if E.title in sl:	
	plotAndSave(E,rng=[0.00,1.00])

	
"""
ROI Analysis
"""
#set measurements to make sure mean, median, area, and position are included
IJ.run("Set Measurements...", "area mean centroid median redirect=None decimal=3");

rm=RoiManager().getInstance()
#choose what central tendancy you want to use Mean or Median
CT='Mean'
if 0!=rm.getCount():
	roiRT=ResultsTable()
	
	rm.runCommand(X1,'Measure')
	tempRT=ResultsTable().getResultsTable()
	roiRT=copyColumn(tempRT,roiRT,'Area','Area')
	roiRT=copyColumn(tempRT,roiRT,'X','X')
	roiRT=copyColumn(tempRT,roiRT,'Y','Y')
	roiRT=copyColumn(tempRT,roiRT,CT,X1.title)
	IJ.run("Clear Results", "");
	
	rm.runCommand(X2,'Measure')
	tempRT=ResultsTable().getResultsTable()
	roiRT=copyColumn(tempRT,roiRT,CT,X2.title)
	IJ.run("Clear Results", "");
	
	rm.runCommand(Y1,'Measure')
	tempRT=ResultsTable().getResultsTable()
	roiRT=copyColumn(tempRT,roiRT,CT,Y1.title)
	IJ.run("Clear Results", "");
	
	rm.runCommand(Y2,'Measure')
	tempRT=ResultsTable().getResultsTable()
	roiRT=copyColumn(tempRT,roiRT,CT,Y2.title)
	IJ.run("Clear Results", "");

	roiFRET(roiRT)
	
	roiRT.show(fnb+' ROI Results')
	roiRT.save(os.path.join(resPath,fnb+' ROI Results.csv'))

	
'''