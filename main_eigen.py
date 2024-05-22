# Importing libraries
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from functions import windowing_individual_breaths, read_csv_files
import os

# fill in patient number (for saving csv filtes)
patient='8_224'

# Loading data 
folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/17_123'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/22_133'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/13_143'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/5_143'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/18_183'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/19_213'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/8_253'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/2_283'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/17_0204'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/7_124'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/16_124'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/2_164'
#folder_path='C:/Users/Jitka/OneDrive - UMCG/Documenten/Stage SEH/Data/8_224'


# read all CSV files in the folder in a list of dataframes 
DF_list=read_csv_files(folder_path)

# extract all capnography and IP data
merged_DF_capno=DF_list[0][(DF_list[0]['Variable']=='Airway CO2 Level Waveform')]
merged_DF_IP=DF_list[0][(DF_list[0]['Variable']=='Impedance Respiration Wave')]

for DF in DF_list:
    DF_capno=DF[(DF['Variable']=='Airway CO2 Level Waveform')]
    DF_IP=DF[(DF['Variable']=='Impedance Respiration Wave')]

    merged_DF_capno=merged_DF_capno.merge(DF_capno, how='outer')
    merged_DF_IP=merged_DF_IP.merge(DF_IP,how='outer')

DF_capno=merged_DF_capno
DF_IP=merged_DF_IP

    
# converting to datetime
DF_capno['Timestamp'] = pd.to_datetime(DF_capno['Timestamp'])
DF_IP['Timestamp'] = pd.to_datetime(DF_IP['Timestamp'])

#ensuring data is in chronological order
DF_capno = DF_capno.sort_values(by='Timestamp')
DF_IP = DF_IP.sort_values(by='Timestamp')

#merge dataframes again to match on time 
merged_for_time=pd.merge_asof(DF_capno,DF_IP,on='Timestamp',by=['Devicetype'],tolerance=pd.Timedelta(seconds=0.16),direction='nearest')

# only keep time instances where capnography and IP data is available
nan_count_df = merged_for_time.isna().sum().sum()
merged_for_time = merged_for_time.dropna()

DF_capno=merged_for_time['Value_x']
DF_IP=merged_for_time['Value_y']

DF_capno=np.array(DF_capno)
DF_IP=np.array(DF_IP)


# creating time axis 

difftime_capno=merged_for_time['Timestamp'].diff().dt.total_seconds()
difftime_capno = difftime_capno.reset_index(drop=True)
difftime_IP=difftime_capno

time_axis_capno=[]

x=difftime_capno[1]

time_axis_capno=np.linspace(0,(len(difftime_capno)-1)*difftime_capno[1],num=(len(difftime_capno)))
time_axis_IP=time_axis_capno


DF_capno=np.char.replace(DF_capno[:].astype(str), ',', '.').astype(float)
DF_IP=np.char.replace(DF_IP[:].astype(str), ',', '.').astype(float)


#plot unfiltered signals
range_tb_displayed=np.arange(int(270/difftime_IP[1]),int(330/difftime_IP[1]))

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(time_axis_capno[:], DF_capno[:])
ax1.set_title('Unfiltered Capnography signal')
ax1.set_xlabel('Time [seconds]')
ax2.plot(time_axis_IP[:], DF_IP[:])
ax2.set_title('Unfiltered Impedance Pneumography signal')
ax2.set_xlabel('Time [seconds]')
plt.show()


## pre processing + count-orig capno

[normalized_capno, normalized_time_capno, peak_values_normalized_capno, peak_timevalues_normalized_capno, troughs_values_normalized_capno, troughs_timevalues_normalized_capno, validpeak_values_capno, validpeak_timevalues_capno,mean_RR_per_window_capno,percentage_high_quality_capno,high_quality_starts_capno,low_quality_starts_capno]=windowing_individual_breaths(DF_capno[:], difftime_capno[1],time_axis_capno,0.25)

# #plotting raw + processed capno
# plt.figure()
# plt.plot(time_axis_capno, DF_capno, 'r', label='Raw signal')
# plt.plot(normalized_time_capno,normalized_capno, 'g', label='Pre-processed signal')
# plt.xlabel('Time [seconds]')
# plt.title('Capnogram')
# plt.legend()
# plt.show()

#plotting peaks and throughs capno + RR intervals
time_axis_RR_per_window_capno=np.arange(0,normalized_time_capno[-1],32)

# fig, (ax1, ax2) = plt.subplots(2, 1)
# ax1.plot(normalized_time_capno,normalized_capno)
# ax1.plot(peak_timevalues_normalized_capno,peak_values_normalized_capno,'ro')
# ax1.plot(troughs_timevalues_normalized_capno,troughs_values_normalized_capno,'bo')
# ax1.plot(validpeak_timevalues_capno,validpeak_values_capno,'ko')
# ax1.set_xlabel('Time [seconds]')
# ax1.set_title('Pre-processed capnography signal and identification of peaks and troughs')
# ax2.step(time_axis_RR_per_window_capno, mean_RR_per_window_capno)
# ax2.set_title('Mean RR intervals per window of 32 seconds')
# ax2.set_xlabel('Time [seconds]')
# ax2.set_ylabel('Time [seconds]')
# plt.show()


## pre processing + count-orig IP
std_threshold=0.25
[normalized_IP, normalized_time_IP, peak_values_normalized_IP, peak_timevalues_normalized_IP, troughs_values_normalized_IP, troughs_timevalues_normalized_IP, validpeak_values_IP, validpeak_timevalues_IP, mean_RR_per_window_IP,percentage_high_quality_IP,high_quality_starts_IP,low_quality_starts_IP]=windowing_individual_breaths(DF_IP[:], difftime_IP[1],time_axis_IP,std_threshold)
    

# # plotting raw + processed IP
# plt.figure()
# plt.plot(time_axis_IP, DF_IP, 'r', label='Raw signal')
# plt.plot(normalized_time_IP,normalized_IP, 'g', label='Pre-processed signal')
# plt.xlabel('Time [seconds]')
# plt.title('Impedance Pneumography')
# plt.legend()
# plt.show()

# #plotting peaks and throughs IP + RR intervals
time_axis_RR_per_window_IP=np.arange(0,normalized_time_IP[-1],32)

# fig, (ax1, ax2) = plt.subplots(2, 1)
# ax1.plot(normalized_time_IP,normalized_IP)
# #ax1.plot(peak_timevalues_normalized_IP,peak_values_normalized_IP,'ro')
# #ax1.plot(troughs_timevalues_normalized_IP,troughs_values_normalized_IP,'bo')
# ax1.plot(validpeak_timevalues_IP,validpeak_values_IP,'ko')
# ax1.set_xlabel('Time [seconds]')
# ax1.set_title('Pre-processed capnography signal and identification of peaks and troughs')
# ax2.step(time_axis_RR_per_window_IP, mean_RR_per_window_IP)
# ax2.set_title('Mean RR intervals per window of 32 seconds')
# ax2.set_xlabel('Time [seconds]')
# ax2.set_ylabel('Time [seconds]')
# plt.show()

# # plotting RR intervals per window for IP and capno

# #difference_mean_RR=abs(mean_RR_per_window_capno-mean_RR_per_window_IP)

# plt.figure()
# plt.step(time_axis_RR_per_window_capno,mean_RR_per_window_capno,where='post',  label='Capnography')
# plt.step(time_axis_RR_per_window_IP, mean_RR_per_window_IP, where='post', label='Impedance Pneumography')
# #plt.step(time_axis_RR_per_window_capno,difference_mean_RR,where='post',label='difference')
# plt.title('Mean RR intervals per window of 32 seconds')
# plt.xlabel('Time [seconds]')
# plt.ylabel('Time [seconds]')
# plt.legend()
# plt.show()





# plotting capno and IP with high and low quality segments

high_quality_starts_capno=np.array(high_quality_starts_capno)
high_quality_starts_IP=np.array(high_quality_starts_IP)

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(normalized_time_capno,normalized_capno, 'r', label='Low quality')
ax2.plot(normalized_time_IP,normalized_IP, 'r',label='Low quality')
ax1.set_title('Pre-processed Capnography')
ax2.set_title('Pre-processed Impedance Pneumography')
ax1.set_xlabel('Time [seconds]')
ax2.set_xlabel('Time [seconds]')


for index,value in enumerate(high_quality_starts_capno):
    high_quality_capno_indices=np.argwhere((normalized_time_capno>=high_quality_starts_capno[index]) & (normalized_time_capno<=high_quality_starts_capno[index]+32))
    ax1.plot(normalized_time_capno[high_quality_capno_indices],normalized_capno[high_quality_capno_indices], 'g', label='High quality')

for index,value in enumerate(high_quality_starts_IP):
    high_qualityp_IP_indices=np.argwhere((normalized_time_IP>=high_quality_starts_IP[index]) & (normalized_time_IP<=high_quality_starts_IP[index]+32))
    ax2.plot(normalized_time_IP[high_qualityp_IP_indices],normalized_IP[high_qualityp_IP_indices],'g',label='High quality')

#ax1.text(10,np.max(normalized_capno)-1,'HQ:'+str(percentage_high_quality_capno) + '%',weight='bold')
#ax2.text(10,np.max(normalized_IP)-1,'HQ:'+str(percentage_high_quality_IP) + '%',weight='bold')

matching_high_quality_starts=[]

for value in high_quality_starts_capno:
    if value in high_quality_starts_IP:
        matching_high_quality_starts.append(value)

matching_high_quality_starts=np.array(matching_high_quality_starts)

matching_HQ=((len(matching_high_quality_starts)*2)/(len(high_quality_starts_capno)+len(low_quality_starts_capno)+len(high_quality_starts_IP)+len(low_quality_starts_IP)))*100
matching_HQ=round(matching_HQ)

print('HQ capno',percentage_high_quality_capno)
print('HQ IP',percentage_high_quality_IP)
print('HQ matched',matching_HQ)

plt.show()

#plotting template/ individual breath cycles + parameter values for a certain segment capno
lower_bound=576      
upperbound=640


range_to_get_specifications=np.argwhere((normalized_time_capno>=lower_bound) & (normalized_time_capno<=upperbound))
specific_valid_peaks_indices=np.argwhere((validpeak_timevalues_capno>=lower_bound) & (validpeak_timevalues_capno<=upperbound))
specific_valid_peaks_indices=specific_valid_peaks_indices.flatten()

valid_peaks_in_range=validpeak_values_capno[specific_valid_peaks_indices]
valid_time_peaks_in_range=validpeak_timevalues_capno[specific_valid_peaks_indices]
normalized_time_capno_range=normalized_time_capno[range_to_get_specifications]
normalized_capno_range=normalized_capno[range_to_get_specifications]




indices_valid_peaks_in_range=np.where(np.isin(normalized_time_capno_range, valid_time_peaks_in_range))[0]

RR_rates=np.diff(valid_time_peaks_in_range)

median_RR_window=np.median(RR_rates)
abnormal_breaths = np.sum((RR_rates > (1.5 * median_RR_window)) | (RR_rates < (0.5 * median_RR_window)))
perc_bad_breaths=abnormal_breaths/(len(RR_rates))

perc_oocupied=(valid_time_peaks_in_range[-1]-valid_time_peaks_in_range[0])/32

mean_RR_window=np.mean(RR_rates)

RR_rate_in_indices=np.diff(indices_valid_peaks_in_range)
RR_rate_in_indices = RR_rate_in_indices[~np.isnan(RR_rate_in_indices)]

mean_RR_in_indices=int(np.mean(RR_rate_in_indices))

std_window=np.std(RR_rates)/mean_RR_window

list_individual_breaths_indices=[]

for peak_index in indices_valid_peaks_in_range:
    lowerlimit_range=(np.round(peak_index-(0.5*mean_RR_in_indices-0.5),0))
    uppelimit_range=(np.round(peak_index+(0.5*mean_RR_in_indices-0.5),0))
    if lowerlimit_range<0 or uppelimit_range>(len(range_to_get_specifications)-1):
            continue
    else:
        peak_range=np.linspace(lowerlimit_range, uppelimit_range, num=(mean_RR_in_indices))
        peak_range=np.array(peak_range)
        list_individual_breaths_indices.append((peak_range))


DF_individual_breaths_indices=pd.DataFrame(list_individual_breaths_indices,index=None, columns=None)

ab_template=[]

for columnindex in DF_individual_breaths_indices.columns:
    column=DF_individual_breaths_indices.iloc[:,columnindex]
    column=np.array(column)
    column=column.astype(int)
    datapoint_template=np.mean(normalized_capno_range[column])
    ab_template.append(datapoint_template)

ab_template=np.array(ab_template)

correlation_coeffients=[]

for index, row in DF_individual_breaths_indices.iterrows():
    validbreath=np.array(row)
    breath_for_row=normalized_capno_range[validbreath.astype(int)]
    breath_for_row=breath_for_row.flatten()
    correlation=np.corrcoef(ab_template,breath_for_row)[0, 1]
    correlation_coeffients.append(correlation)

mean_corr_coeff=np.mean(np.array(correlation_coeffients))

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(normalized_time_capno[range_to_get_specifications],normalized_capno[range_to_get_specifications], 'red')
ax1.plot(valid_time_peaks_in_range,valid_peaks_in_range,'ko')
ax1.set_title('Low quality Capnography segment')
ax1.set_xlabel('Time [seconds]', fontsize=13)



for index, row in DF_individual_breaths_indices.iterrows():
     validbreath=np.array(row)
     time_segment=np.linspace(0,np.size(validbreath)*0.2, num=np.size(validbreath))
     ax2.plot(time_segment,normalized_capno_range[validbreath.astype(int)],'black')

ax2.plot(time_segment,ab_template, 'red')

ax2.set_title('Template')
ax2.set_xlabel('Time [seconds]')

plt.show()

print('RR-intervals',RR_rates)
print('Percentage occupied', perc_oocupied)
print('Percentage bad breaths', perc_bad_breaths)
print('Standard deviations RR-intervals', std_window)
print('Mean correlation coefficient',mean_corr_coeff )

#plotting template/ individual breath cycles + parameter values for a certain segment IP
lower_bound=576      
upperbound=640

range_to_get_specifications=np.argwhere((normalized_time_IP>=lower_bound) & (normalized_time_IP<=upperbound))
specific_valid_peaks_indices=np.argwhere((validpeak_timevalues_IP>=lower_bound) & (validpeak_timevalues_IP<=upperbound))
specific_valid_peaks_indices=specific_valid_peaks_indices.flatten()

valid_peaks_in_range=validpeak_values_IP[specific_valid_peaks_indices]
valid_time_peaks_in_range=validpeak_timevalues_IP[specific_valid_peaks_indices]
normalized_time_IP_range=normalized_time_IP[range_to_get_specifications]
normalized_IP_range=normalized_IP[range_to_get_specifications]


indices_valid_peaks_in_range=np.where(np.isin(normalized_time_IP_range, valid_time_peaks_in_range))[0]

RR_rates=np.diff(valid_time_peaks_in_range)

median_RR_window=np.median(RR_rates)
abnormal_breaths = np.sum((RR_rates > (1.5 * median_RR_window)) | (RR_rates < (0.5 * median_RR_window)))
perc_bad_breaths=abnormal_breaths/(len(RR_rates))

perc_oocupied=(valid_time_peaks_in_range[-1]-valid_time_peaks_in_range[0])/32

mean_RR_window=np.mean(RR_rates)

RR_rate_in_indices=np.diff(indices_valid_peaks_in_range)
RR_rate_in_indices = RR_rate_in_indices[~np.isnan(RR_rate_in_indices)]

mean_RR_in_indices=int(np.mean(RR_rate_in_indices))

std_window=np.std(RR_rates)/mean_RR_window

list_individual_breaths_indices=[]

for peak_index in indices_valid_peaks_in_range:
    lowerlimit_range=(np.round(peak_index-(0.5*mean_RR_in_indices-0.5),0))
    uppelimit_range=(np.round(peak_index+(0.5*mean_RR_in_indices-0.5),0))
    if lowerlimit_range<0 or uppelimit_range>(len(range_to_get_specifications)-1):
            continue
    else:
        peak_range=np.linspace(lowerlimit_range, uppelimit_range, num=(mean_RR_in_indices))
        peak_range=np.array(peak_range)
        list_individual_breaths_indices.append((peak_range))


DF_individual_breaths_indices=pd.DataFrame(list_individual_breaths_indices,index=None, columns=None)

ab_template=[]

for columnindex in DF_individual_breaths_indices.columns:
    column=DF_individual_breaths_indices.iloc[:,columnindex]
    column=np.array(column)
    column=column.astype(int)
    datapoint_template=np.mean(normalized_IP_range[column])
    ab_template.append(datapoint_template)

ab_template=np.array(ab_template)

correlation_coeffients=[]

for index, row in DF_individual_breaths_indices.iterrows():
    validbreath=np.array(row)
    breath_for_row=normalized_IP_range[validbreath.astype(int)]
    breath_for_row=breath_for_row.flatten()
    correlation=np.corrcoef(ab_template,breath_for_row)[0, 1]
    correlation_coeffients.append(correlation)

mean_corr_coeff=np.mean(np.array(correlation_coeffients))

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(normalized_time_IP[range_to_get_specifications],normalized_IP[range_to_get_specifications] , 'green')
ax1.plot(valid_time_peaks_in_range,valid_peaks_in_range,'ko')
ax1.set_title('LQ Impedance Pneumography segment', fontsize=16 )
ax1.set_xlabel('Time [seconds]', fontsize=13)



for index, row in DF_individual_breaths_indices.iterrows():
     validbreath=np.array(row)
     time_segment=np.linspace(0,np.size(validbreath)*0.2, num=np.size(validbreath))
     ax2.plot(time_segment, normalized_IP_range[validbreath.astype(int)],'black')

ax2.plot(time_segment,ab_template, 'red')

ax2.set_title('Template')
ax2.set_xlabel('Time (in s)')

plt.show()

print('RR-intervals',RR_rates)
print('Percentage occupied', perc_oocupied)
print('Percentage bad breaths', perc_bad_breaths)
print('Standard deviations RR-intervals', std_window)
print('Mean correlation coefficient',mean_corr_coeff )

# extract starting points 
high_quality_starts_capno_diff=high_quality_starts_capno/32
high_quality_starts_capno_diff = high_quality_starts_capno_diff.astype(int)

matching_high_quality_starts_diff=matching_high_quality_starts/32
matching_high_quality_starts_diff = matching_high_quality_starts_diff.astype(int)

# convert RR interval arrays to dataframes for high capno data 
DF_RR_intervals_capno_HQ_capno = pd.DataFrame(mean_RR_per_window_capno[high_quality_starts_capno_diff])
DF_RRintervals_IP_HQ_capno= pd.DataFrame(mean_RR_per_window_IP[high_quality_starts_capno_diff])

# convert RR interval arrays to dataframes for matching high quality
DF_RR_intervals_capno_matching = pd.DataFrame(mean_RR_per_window_capno[matching_high_quality_starts_diff])
DF_RR_intervals_IP_matching= pd.DataFrame(mean_RR_per_window_IP[matching_high_quality_starts_diff])
  

# create dataframe for capnography indicating its quality (only capno and matching)
qualitycolumn=[]

for i in np.arange(0, len(mean_RR_per_window_capno)):
    if i in matching_high_quality_starts_diff:
        qualitycolumn.append('HQ matching')
    elif i in high_quality_starts_capno_diff:
        qualitycolumn.append('HQ capnography')
    else:
        qualitycolumn.append('Low quality')

mean_RR_per_window_capno = mean_RR_per_window_capno.tolist()


capnography=pd.DataFrame({'RR interval':mean_RR_per_window_capno,'Quality':qualitycolumn })

# create dataframe for IP indicating its quality (only capno and matching)



IP=pd.DataFrame({'RR interval':mean_RR_per_window_IP,'Quality':qualitycolumn })

# save the dataframes containing RR intervals in csv 
directory_path='C:/Users/Jitka/OneDrive/Documenten/Technische Geneeskunde Master/Stages/Stage 2  - SEH UMCG/Scripts/'
folder_name=str(std_threshold)
folder_path_first = os.path.join(directory_path, folder_name)

if not os.path.exists(folder_path_first):
    os.makedirs(folder_path_first)

folder_path = os.path.join(folder_path_first, 'capnography')

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

filename = "RR_intervals_capno_" + patient + ".csv"
full_path = os.path.join(folder_path, filename)
capnography.to_csv(full_path)


folder_path = os.path.join(folder_path_first, 'IP')

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

filename = "RR_intervals_IP_" + patient + ".csv"
full_path = os.path.join(folder_path, filename)
# IP.to_csv(full_path)



