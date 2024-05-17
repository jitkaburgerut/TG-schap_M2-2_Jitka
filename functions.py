import numpy as np
from scipy.signal.windows import tukey
from scipy.signal import butter,filtfilt, detrend
from scipy.interpolate import interpolate
import os
import pandas as pd
import math
import matplotlib.pyplot as plt

def read_csv_files(folder_path):
    DFs = []
    
    for file in os.listdir(folder_path):
        if file.endswith(".csv"):
            file_path = os.path.join(folder_path, file)
            DFs.append(pd.read_csv(file_path, sep=';', decimal=','))
    return DFs

def read_csv_files_BA(folder_path):
    DFs = []
    
    for file in os.listdir(folder_path):
        if file.endswith(".csv"):
            file_path = os.path.join(folder_path, file)
            DFs.append(pd.read_csv(file_path))
    return DFs

def count_orig(signal):

    #local extrema
    # wanneer platte pieken --> niet meegenomen

    localpeaksindices=[]
    localtroughsindices=[]
    
    for count, value in enumerate(signal):
        if count==0:
            if signal[count]>signal[count+1] or signal[count]==max(signal):
                localpeaksindices.append(count)
            elif signal[count]<signal[count+1] or min(signal):
                localtroughsindices.append(count)
        elif (count==len(signal)-1):
            if signal[count]>signal[count-1] or signal[count]==max(signal):
                localpeaksindices.append(count)
            elif signal[count]<signal[count-1] or min(signal):
                localtroughsindices.append(count)
        else:
            if (signal[count]>signal[count-1] and signal[count]>signal[count+1]) or signal[count]==max(signal) :
                localpeaksindices.append(count)
            elif (signal[count]<signal[count-1] and signal[count]<signal[count+1]) or signal[count]==min(signal):
                localtroughsindices.append(count)
    
    localpeaks=signal[localpeaksindices]
    localtroughs=signal[localtroughsindices]

    threshold_peaks=np.quantile(localpeaks,0.75)*0.2
    threshold_troughs=np.quantile(localtroughs,0.25)*0.2

    relevantpeakindices=[]
    relevanttroughindices=[]

    for count, value in enumerate(localpeaks):
        if value>threshold_peaks:
            relevantpeakindices.append(localpeaksindices[count])
        
    for count, value in enumerate(localtroughs):
        if value<threshold_troughs:
            relevanttroughindices.append(localtroughsindices[count])

    #valid breaths
            
    valid_peak_indices=[]
    lowerbound=0  

    for count0, value0 in enumerate(relevanttroughindices):
            

            
            if count==(len(relevanttroughindices)-1):
                peaks_indices_last_breath=[value for value in relevantpeakindices if lowerbound < value < value0]
                if len(peaks_indices_last_breath)>0:
                    peak_index_last_breath=max(peaks_indices_last_breath, key=lambda i: signal[i])
                    valid_peak_indices.append(peak_index_last_breath)
                lowerbound=value0 
                # peaks_indices_breath=[value for value in relevantpeakindices if lowerbound < value < 20000]
                # if len(peak_index_breath)>0:
                #     peak_index_breath=max(peaks_indices_breath, key=lambda i: signal[i])
                #     valid_peak_indices.append(peak_index_breath)
           
            else: # inbouwen dat er ook geen peak indices meer kunnen zijn 
                peaks_indices_breath=[value for value in relevantpeakindices if lowerbound < value < value0]
                if len(peaks_indices_breath)>0:
                    peak_index_breath=max(peaks_indices_breath, key=lambda i: signal[i])
                    valid_peak_indices.append(peak_index_breath) #klopt niet 
                lowerbound=value0


    return relevantpeakindices, relevanttroughindices,valid_peak_indices




def butter_lowpass_filter(data, difftime): # change to kaiserord?
    #data=data.fillna(method='ffill')
    alpha=4/(len(data)*difftime)
    tukey_window=tukey(len(data),alpha=alpha)
    tukey_signal=detrend((data)*tukey_window)

    order=2
    cutoff=1
    fs=1/(difftime)
    nyq = 0.5 * fs 
    normal_cutoff = cutoff / nyq
   
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, tukey_signal)
    return y

def windowing_individual_breaths(data, difftime, timeaxis, std):
    
    lowpass_filtered=butter_lowpass_filter(data, difftime)

    windowlength=round(32/difftime)
    window_start_indices=np.arange(0,len(lowpass_filtered),windowlength)

    window_starts=np.arange(0,timeaxis[-1],32)
    window_starts = window_starts.astype(int)

    mean_RR_per_window=[]
    norm_std_RR=[]
    normalized=[]
    normalized_time=[]
    peak_values_normalized=[]
    troughs_values_normalized=[]
    peak_timevalues_normalized=[]
    troughs_timevalues_normalized=[]
    validpeak_values=[]
    validpeak_timevalues=[]

    high_quality_starts=[]
    low_quality_starts=[]

    len_window_starts=len(window_starts)

    STDi=[]
    breath_duration=[]
    occupation_vb=[]
    mean_correlation=[]


    for window_number in range(len_window_starts-1): #laatste window wordt dus niet meegenomen
        win_start=int(window_starts[window_number])
        win_end=int(win_start+32)
        if window_number==0:
            window_indices=np.argwhere((timeaxis>=win_start) & (timeaxis<=win_end))
        else:
            window_indices=np.argwhere((timeaxis>win_start) & (timeaxis<=win_end))
        window_indices=window_indices.flatten()
        window=lowpass_filtered[window_indices]
        window_time=timeaxis[window_indices]

        #lineaire interpoleratie voor downsampling naar 5 Hz
        interp_data_t=[]

        if window_number==0:
            sample_numbers=int((32/0.2)+1)
            interp_data_t= np.linspace(win_start,win_end,num=sample_numbers)

        else: 
            sample_numbers=int(32/0.2)
            interp_data_t= np.linspace((win_start+0.2),win_end,num=sample_numbers)

        interp_func = interpolate.interp1d(window_time, window, kind='linear')
        inter_pol = interp_func(interp_data_t)

        std_window=np.std(inter_pol)
        mean_window=np.mean(inter_pol)
        normalized_window=(inter_pol-mean_window)/std_window

        normalized.extend(normalized_window)
        normalized_time.extend(interp_data_t)

        [peaksindices_window,troughsindices_window,valid_peak_indices]=count_orig(normalized_window[:])

        peak_values_normalized.extend(normalized_window[peaksindices_window])
        peak_timevalues_normalized.extend(interp_data_t[peaksindices_window])
        troughs_values_normalized.extend(normalized_window[troughsindices_window])
        troughs_timevalues_normalized.extend(interp_data_t[troughsindices_window])
        validpeak_values.extend(normalized_window[valid_peak_indices])
        validpeak_timevalues.extend(interp_data_t[valid_peak_indices])

        #SQI
        valid_peak_tvalues=interp_data_t[valid_peak_indices]
        valid_peak_tvalues=np.array(valid_peak_tvalues)

        RR_rates=np.diff(valid_peak_tvalues)
        mean_RR_window=np.mean(RR_rates)
        

        RR_rate_in_indices=np.diff(valid_peak_indices)

        RR_rate_in_indices = RR_rate_in_indices[~np.isnan(RR_rate_in_indices)]

        # if np.isnan(RR_rate_in_indices).any():
        #     break

        # if not len(RR_rate_in_indices):
        #     mean_RR_per_window.append(np.nan)
    
        #     continue

        if len(RR_rate_in_indices)==0:
            mean_RR_per_window.append(np.nan)
            continue

        mean_RR_in_indices=int(np.mean(RR_rate_in_indices))
        
       
       #normalized std
        norm_std_RR_window=np.std(RR_rates)/mean_RR_window #coefficient of variation 
        
        

        

        #abnormal breath duration
        median_RR_window=np.median(RR_rates)
        abnormal_breaths = np.sum((RR_rates > (1.5 * median_RR_window)) | (RR_rates < (0.5 * median_RR_window)))
        condition=(RR_rates > (1.5 * median_RR_window)) | (RR_rates < (0.5 * median_RR_window))
        indices_bad_breaths=np.where(condition)[0]
        indices_bad_breaths = np.array(indices_bad_breaths, dtype=np.int32)
        perc_bad_breaths=abnormal_breaths/(len(RR_rates))

        perc_oocupied=(valid_peak_tvalues[-1]-valid_peak_tvalues[0])/32

        #extra regel 

        # if perc_bad_breaths<=0.20:
            
        #     RR_rates = np.delete(RR_rates, indices_bad_breaths)
        #     valid_peak_indices = np.delete(valid_peak_indices, indices_bad_breaths)


        #     mean_RR_window=np.mean(RR_rates)

        #     RR_rate_in_indices=np.diff(valid_peak_indices)

        #     RR_rate_in_indices = RR_rate_in_indices[~np.isnan(RR_rate_in_indices)]

        #     mean_RR_in_indices=int(np.mean(RR_rate_in_indices))

        #     norm_std_RR_window=np.std(RR_rates)/mean_RR_window

        #     abnormal_breaths = np.sum((RR_rates > (1.5 * median_RR_window)) | (RR_rates < (0.5 * median_RR_window)))
        #     perc_bad_breaths=abnormal_breaths/(len(RR_rates))

        #     perc_oocupied=(np.sum(RR_rates))/32

        #     median_RR_window=np.median(RR_rates)

            


        mean_RR_per_window.append(mean_RR_window)
        norm_std_RR.append(norm_std_RR_window)    
            


        

    
        

        #creating template

        list_individual_breaths_indices=[]

        for peak_index in valid_peak_indices:
            lowerlimit_range=(np.round(peak_index-(0.5*mean_RR_in_indices-0.5),0))
            uppelimit_range=(np.round(peak_index+(0.5*mean_RR_in_indices-0.5),0))
            if lowerlimit_range<0 or uppelimit_range>(len(normalized_window)-1):
                continue
            else:
                peak_range=np.linspace(lowerlimit_range, uppelimit_range, num=(mean_RR_in_indices))
                peak_range=np.array(peak_range)
                list_individual_breaths_indices.append((peak_range))




        DF_individual_breaths_indices=pd.DataFrame(list_individual_breaths_indices,index=None, columns=None)

        # plt.figure()

        # #plt.plot(normalized_window)
        # #plt.plot(valid_peak_indices, normalized_window[valid_peak_indices], 'go')

        # for index, row in DF_individual_breaths_indices.iterrows():
        #     validbreath=np.array(row)
        #     plt.plot(normalized_window[validbreath.astype(int)])

        # plt.show()

    
        ab_template=[]

        for columnindex in DF_individual_breaths_indices.columns:
            column=DF_individual_breaths_indices.iloc[:,columnindex]
            column=np.array(column)
            column=column.astype(int)
            datapoint_template=np.mean(normalized_window[column])
            ab_template.append(datapoint_template)

        ab_template=np.array(ab_template)

        correlation_coeffients=[]

        for index, row in DF_individual_breaths_indices.iterrows():
            validbreath=np.array(row)
            breath_for_row=normalized_window[validbreath.astype(int)]
            correlation=np.corrcoef(ab_template,breath_for_row)[0, 1]
            correlation_coeffients.append(correlation)

        mean_corr_coeff=np.mean(np.array(correlation_coeffients))

        # print('Breath durations',RR_rates)
        # print('Std', norm_std_RR_window)
        # print('Mean correlation', mean_corr_coeff)
        # print('total breath duration + occupied', np.sum(RR_rates),perc_oocupied)
        # print('num bad breaths+perc bad breaths',abnormal_breaths, perc_bad_breaths)

        # fig, (ax1, ax2) = plt.subplots(2, 1)
        # ax1.plot(ab_template, 'red')

        # for index, row in DF_individual_breaths_indices.iterrows():
        #     validbreath=np.array(row)
        #     ax1.plot(normalized_window[validbreath.astype(int)],'black')
        

        # ax2.plot(interp_data_t,normalized_window)
        # ax2.plot(interp_data_t[valid_peak_indices],normalized_window[valid_peak_indices],'ko')

        plt.show()

        # determine high/low quality

        
        STDi.append(norm_std_RR_window)
        breath_duration.append(perc_bad_breaths)
        occupation_vb.append(perc_oocupied)
        mean_correlation.append(mean_corr_coeff)

        if perc_oocupied>0.6 and perc_bad_breaths<0.15 and norm_std_RR_window<std and mean_corr_coeff>0.75:
            high_quality_starts.append(int(window_starts[window_number]))
        else:
            low_quality_starts.append(int(window_starts[window_number]))

        #std deviation has the most influence

    l_STD= [i for i in STDi if i>0.25]
    num_STD=len(l_STD) 

    l_breath_duration= [i for i in breath_duration if i>0.15]
    num_breath_duration=len(l_breath_duration)

    l_occupation_vb=[i for i in occupation_vb if i<0.60]
    num_o_vb=len(l_occupation_vb)

    l_mean_cor=[i for i in mean_correlation if i<0.75]
    num_mean_correlation=len(l_mean_cor)

    normalized=np.array(normalized)
    normalized_time=np.array(normalized_time)
    peak_values_normalized=np.array(peak_values_normalized)
    peak_timevalues_normalized=np.array(peak_timevalues_normalized)
    troughs_values_normalized=np.array(troughs_values_normalized)
    troughs_timevalues_normalized=np.array(troughs_timevalues_normalized)
    validpeak_values=np.array(validpeak_values)
    validpeak_timevalues=np.array(validpeak_timevalues)
    mean_RR_per_window=np.array(mean_RR_per_window)

    percentage_high_quality=round((len(high_quality_starts)/(len(high_quality_starts)+len(low_quality_starts)))*100)

    print('number HQ', len(high_quality_starts))
    print('number LQ', len(low_quality_starts))
    print('number std', num_STD)
    print('number breath duration', num_breath_duration)
    print('number occupied valid breaths', num_o_vb)
    print('number mean correlation', num_mean_correlation)

    return normalized, normalized_time, peak_values_normalized, peak_timevalues_normalized, troughs_values_normalized, troughs_timevalues_normalized, validpeak_values, validpeak_timevalues,mean_RR_per_window, percentage_high_quality, high_quality_starts,low_quality_starts

def BlandAltman(data1,data2):
    x_axis=[]
    y_axis=[]

    for i in range(0,len(data1),1):
        mean=np.mean([data1[i],data2[i]])
        x_axis.append(mean)
        difference=(data2[i]-data1[i])
        y_axis.append(difference) 

    x_axis=np.array(x_axis)
    y_axis=np.array(y_axis)
    y_axis = y_axis[~np.isnan(y_axis)]
    x_axis = x_axis[~np.isnan(x_axis)]

    mean=np.mean(y_axis)
    
    std_diff = np.std(y_axis)  # ddof=1 for sample standard deviation

    # Calculate limits of agreement
    limit_of_agreement_plus = 1.96 * std_diff +mean
    limit_of_agreement_min =-1.96 * std_diff+ mean



    return x_axis, y_axis, limit_of_agreement_min,limit_of_agreement_plus