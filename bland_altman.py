import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from functions import read_csv_files_BA
import pingouin as pg

folder_path_capno='C:/Users/Jitka/OneDrive/Documenten/Technische Geneeskunde Master/Stages/Stage 2  - SEH UMCG/Scripts/0.25/capnography'
folder_path_IP='C:/Users/Jitka/OneDrive/Documenten/Technische Geneeskunde Master/Stages/Stage 2  - SEH UMCG/Scripts/0.25/IP'
DF_list_capno=read_csv_files_BA(folder_path_capno)
DF_list_IP=read_csv_files_BA(folder_path_IP)

all_capno=0
capno_HQ_matching=[]
capno_HQ_capno=[]

for DF in DF_list_capno:
    capno_all_patient=len(DF)
    RR_capno_matchingHQ_patient= DF.loc[DF['Quality'] == 'HQ matching', 'RR interval']
    RR_capno_capnoHQ_patient= DF.loc[DF['Quality'] == 'HQ capnography', 'RR interval']
    
    all_capno=all_capno+capno_all_patient
    capno_HQ_matching.extend(RR_capno_matchingHQ_patient)
    capno_HQ_capno.extend(RR_capno_capnoHQ_patient)

all_IP=0
IP_HQ_matching=[]
IP_HQ_capno=[]

for DF in DF_list_IP:
    IP_all_patient=len(DF)
    RR_IP_matchingHQ_patient= DF.loc[DF['Quality'] == 'HQ matching', 'RR interval']
    RR_IP_capnoHQ_patient= DF.loc[DF['Quality'] == 'HQ capnography', 'RR interval']
    
    all_IP=all_IP+IP_all_patient
    IP_HQ_matching.extend(RR_IP_matchingHQ_patient)
    IP_HQ_capno.extend(RR_IP_capnoHQ_patient)




capno_HQ_matching=np.array(capno_HQ_matching)
capno_HQ_capno=np.array(capno_HQ_capno)
IP_HQ_matching=np.array(IP_HQ_matching)
IP_HQ_capno=np.array(IP_HQ_capno)

meanRR_capno=60/np.mean(capno_HQ_capno)
meanRR_IP=60/np.mean(IP_HQ_capno)

stdRR_capno=np.std(60/capno_HQ_capno)
stdRR_IP=np.std(60/IP_HQ_capno)

capno_all_HQ_capno=np.concatenate((capno_HQ_matching,capno_HQ_capno))
IP_all_HQ_capno=np.concatenate((IP_HQ_matching,IP_HQ_capno))

x_axis=[]
y_axis=[]

def BlandAltman(data1,data2):
    
    for i in range(0,len(data1),1):
        mean=np.mean([60/data1[i],60/data2[i]])
        x_axis.append(mean)
        difference=(60/data2[i]-60/data1[i])
        y_axis.append(difference) 

    mean=np.mean(y_axis)
    
    std_diff = np.std(y_axis)  # ddof=1 for sample standard deviation

    # Calculate limits of agreement
    limit_of_agreement_plus = 1.96 * std_diff +mean
    limit_of_agreement_min =-1.96 * std_diff+ mean

    return limit_of_agreement_min,limit_of_agreement_plus
        

         
[limitmin,limitplus]=BlandAltman(capno_HQ_capno,IP_HQ_capno)

linemin=np.ones(len(y_axis))*limitmin
lineplus=np.ones(len(y_axis))*limitplus
meanline=np.ones(len(y_axis))*np.mean(y_axis)

x_axis=np.array(x_axis)
y_axis=np.array(y_axis)


#plotting
bin_width=1
data_range = np.max(y_axis) - np.min(y_axis)
num_bins = int(data_range / bin_width)


fig, (ax1, ax2,ax3) = plt.subplots(1, 3)
ax1.scatter((60/np.array(capno_HQ_capno)),(60/np.array(IP_HQ_capno)))
ax1.plot([0,35],[0,35], 'k')
ax1.set_ylabel('IP measurement (in bpm)', fontsize=13)
ax1.set_xlabel('Capnography measurement (in bpm)', fontsize=13)

ax2.scatter(x_axis,y_axis, color='blue')
ax2.plot(x_axis,linemin, 'red')
ax2.plot(x_axis,lineplus,'red')
ax2.plot(x_axis,meanline    ,'black')
ax2.set_ylabel('Difference (in bpm)', fontsize=13)
ax2.set_xlabel('Mean (in bpm)', fontsize=13)

ax3.hist((y_axis),bins=num_bins)
#ax3.title('Differences between capnography and IP for HQ capnography data')
ax3.set_xlabel('Difference (in bpm)', fontsize=13)
ax3.set_ylabel('Count', fontsize=13)



fig.suptitle('Data Analysis LQ IP data',fontsize=16)
plt.show()







targets1=np.arange(1,(len(capno_HQ_capno)+1))

targets=np.concatenate((targets1, targets1))

Ratings=np.concatenate((capno_HQ_capno,IP_HQ_capno))

plt.show()
# Create the list
Y = ['Y'] * len(targets1)
X = ['X'] * len(targets1)
Raters = Y + X

df = pd.DataFrame({'Targets': targets,
                   'Rater': Raters,
                   'Ratings': Ratings})

icc = pg.intraclass_corr(data=df,targets='Targets',raters='Rater',  ratings='Ratings').round(3)
icc.set_index("Type")

# sensitivity and specificity

TP=0
FP=0
TN=0
FN=0

for index,value in enumerate(capno_HQ_matching):
    RR_freq_capno_matching=60/capno_HQ_matching[index]
    RR_freq_IP_matching=60/IP_HQ_matching[index]
    difference_matching=RR_freq_capno_matching-RR_freq_IP_matching
    if abs(difference_matching)<=0.5:
        TP=TP+1
    else:
        FP=FP+1
    
for index,value in enumerate(capno_HQ_capno):
    RR_freq_capno_capno=60/capno_HQ_capno[index]
    RR_freq_IP_capno=60/IP_HQ_capno[index]
    difference_capno=RR_freq_capno_capno-RR_freq_IP_capno
    if abs(difference_capno)<=0.5:
        FN=FN+1
    else:
        TN=TN+1



sensitiviteit=TP/(TP+FN)
specificiteit=TN/(TN+FP)
preciscion=TP/(TP+FP)
NPV=TN/(TN+FN)
accuracy=(TP+TN)/(TN+TP+FP+FN)



Matched_perc=(len(capno_HQ_matching*2))/(all_IP+all_capno)
HQ_perc_capno=(len(capno_HQ_matching)+len(capno_HQ_capno))/(all_capno)
                                                         
x=1