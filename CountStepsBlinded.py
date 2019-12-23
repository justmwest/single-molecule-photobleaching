#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 11:12:11 2019

@author: Justin
"""

#Outline for this script
#Step 1: Find this notebook
#Step 2: Make assignment folders in all subdirectories (1 step, etc.)
#Step 3: Make a matrix with the dimensions of the directory (expts, movies)
#Step 4: Make a randomly ordered list of the coordinates of each movie
#Step 5: Start list, open graph, prompt user to input count
#Step 6: Copy the data to the corresponding subdirectory, close graph
#Step 7: Repeat steps 5 & 6 until end of list
#Step 8: Show a table of the final results

import glob
import os
import random
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from tabulate import tabulate


######### Step 1 #########
### Find this notebook ###
##########################





#set path for the day, not the slide/experiment
path = "/Users/Justin/Documents/Projects/InProgress/TM-EphA2-PIP2/Data/SMALPs/SM bleaching/2019-12-20 copy/"
os.chdir(path);
dir_names = sorted(next(os.walk(path))[1])
n_dir = len(dir_names) #get the number of subdirectories in path
print(n_dir," folders")





#gen a list of the number of movies in each subdirectory
extension = 'traces'
def get_n_movies(directory):
    '''
    input:
        int, directory index starting from 0
    output:
        returns the number of movies in the directory
    '''
    os.chdir(dir_names[directory])
    tmp = len([m for m in glob.glob('*.{}'.format(extension))])
    os.chdir(path)
    return tmp   

n_movies=[]    
for i in range(len(dir_names)):
    n_movies.append(get_n_movies(i))
print(n_movies," movies")

sum_movies = 0
for i in n_movies:
    sum_movies += i
print(sum_movies," movies total")





def get_n_traces(directory,movie):
    '''
    input:
        int, directory index starting from 0
        int, movie index starting from 1
    output:
        returns the number of traces in a movie
    '''
    os.chdir(path)
    os.chdir(dir_names[directory])
    fname = "hel"+str(movie)+".traces"
    fid = open(fname,'r') # the open file object
    length = int(np.fromfile(fid,'int32', 1)) # number of timesteps, must keep.
    n_traces = int(np.fromfile(fid, 'int16', 1)) # number of traces
    n_traces = int(n_traces/2)
    fid.close()
    os.chdir(path)
    return n_traces

grand_total_traces = 0
for j in range(len(dir_names)):
    for i in range(1,n_movies[j]+1):
        n_traces = get_n_traces(j,i)
        grand_total_traces += n_traces
print("The total number of molecules collected is "+str(grand_total_traces))






def get_coordinates(number):
    '''
    Returns coordinates for the molecule number entered by
    making a list of coordinates of len grand_total_traces
    input:
        int, between 0 and grand_total_traces
    output:
        list, coordinates [directory, movie, trace]
    Dependencies:
        function, get_n_traces
    Notes:
        -Perhaps could be condensed using unravel
    '''
    top_lvl_lofl = []
    for j in range(n_dir):
        dir_lvl_lofl = []
        for i in range(1,n_movies[j]+1):
            list1 = [j]*get_n_traces(j,i) # directory number
            list2 = [i]*get_n_traces(j,i) # movie number
            list3 = list(range(1,get_n_traces(j,i)+1)) # trace number
            #put them all together
            lofl = []
            for i in range(get_n_traces(j,i)):
                lofl.append([list1[i],list2[i],list3[i]])
            dir_lvl_lofl.append(lofl)
        dir_lvl_lofl = [item for sublist in dir_lvl_lofl for item in sublist] #Flatten one level
        top_lvl_lofl.append(dir_lvl_lofl)
        #print(dir_lvl_lofl)
    top_lvl_lofl = [item for sublist in top_lvl_lofl for item in sublist] #Flatten one level
    #print(len(top_lvl_lofl)) #test, should match n of molecules
    return top_lvl_lofl[number]


    
def gen_master_list():
    '''
    Makes a counting list of len grand_total_traces,
    Then shuffles it to make the order of files
    '''
    tmp_list = list(range(grand_total_traces))
    random.shuffle(tmp_list)
    return tmp_list
order = gen_master_list()





##########################
#### Make Directories ####
##########################
#   I would like to make a new directory tree for each analysis.
#   That way, I don't have to worry about data loss when starting
#   a new analysis. 


# Make parent analyses folder if not already present
os.chdir(path)
directory_basename = os.path.basename(os.path.normpath(path)) #returns folder name
analyses_dirname = directory_basename+" analyses"
os.chdir(os.pardir)
if not os.path.exists(analyses_dirname):
    os.mkdir(analyses_dirname)
os.chdir(analyses_dirname)
analyses_path = os.getcwd()


# Make a "current analysis" folder
now = datetime.now() #current time for analysis folder
current_time = now.strftime("%Y-%m-%d_%H-%M-%S")
print("Current Time =", current_time)
os.chdir(analyses_path)
os.mkdir(current_time)
os.chdir(current_time)
current_analysis_path=os.getcwd()

# Make directory tree for analysis files        
for j in range(len(dir_names)):
    for i in range(1,n_movies[j]+1):
        os.makedirs(dir_names[j]+"/hel"+str(i)+"/Rejected")
        for step in range(1,9):
            if step == 1:
                suffix = " step"
            else:
                suffix = " steps"
            os.makedirs(dir_names[j]+"/hel"+str(i)+"/"+str(step)+suffix)



##########################
####### Read data ########
##########################
#   Here is where I make the functions to use in the main analysis.
            
timeunit=0.1
ymin = -50
ymax = 1000
leakage=0


def read_data(coordinates):
    '''
    Master file reading function
    Assigns global variables that can be used in other functions
    
    input:
        coordinates as list. [directory, movie, trace]
    output:
        null
    '''
    global directory, movie, trace
    directory = coordinates[0]
    movie = coordinates[1]
    trace = coordinates[2]

    global fname, fid, length, n_traces, raw
    os.chdir(dir_names[directory])
    fname = "hel"+str(movie)+".traces"
    fid = open(fname,'r') # the open file object
    length = int(np.fromfile(fid,'int32', 1)) # number of timesteps
    n_traces = int(np.fromfile(fid, 'int16', 1)) # number of traces * 2
    raw = np.fromfile(fid, 'int16', n_traces*length);
    fid.close()
    
    global index, data
    index = range(n_traces*length) # list of indeces
    data = np.zeros((n_traces, length)) #Data=zeros(Ntraces,len);
    data[np.unravel_index(index, data.shape, 'F')] = raw[index] #StkOvflw usr "Suever"

    global donor, acceptor, time 
    donor = np.zeros((n_traces//2, length)) #gen empty array for donor
    acceptor = np.zeros((n_traces//2, length)) #gen empty array for acceptor
    for i in range(0,int(n_traces),2): #populate donor and acceptor arrays
       donor[i//2] = data[i]
       acceptor[i//2] = data[i+1]
    time = list(range(length))    
    time = [x * timeunit for x in time] # make a time list for the x values
    os.chdir(path)
    return



def plot(blinded=True):
    '''
    Plotting function
    
    input:
        null
    output:
        a plot using global vars assigned in read_data function.
    '''
    plt.figure(figsize=(9,3))
    plt.plot(time, acceptor[trace-1]-leakage*donor[trace-1], 'r')
    plt.plot(time, donor[trace-1], 'g')
    axes = plt.gca() #stands for get current axes. from stkovflw usr Hima.
    axes.set_xlim([0,max(time)]) # set x axis range
    axes.set_ylim([ymin,ymax]) # set y axis range
    plt.grid(True)
    plt.ylabel('Intensity')
    plt.xlabel('Time (seconds, frames x '+str(timeunit)+')')
    if not blinded:
        plt.title('Folder: '+str(dir_names[directory])+'. Movie: '+str(movie)+'. Molecule: '+str(trace)) # unblinded only
    plt.show()


 
def export(data_list,file_name="untitled.xlsx"):
    '''
    Exports data to excel file in state 1 2 3 format.
    Works even if a value is none!
    
    input:
        data_list, list of lists of ints or floats
        file_name, str. include the extension
    output:
        saves an excel file in the cwd
    '''
    list_for_export=[]
    for each in data_list:
        if each is None:
            list_for_export.append(np.ndarray([]).flatten())
        else:
            list_for_export.append(each.flatten())
    df = pd.DataFrame(list_for_export, index=None)
    df = df.T
    df.to_excel(file_name,index=None,header=None)
    
    
def has_numbers(inputString):
    '''
    returns True if the string contains numbers
    '''
    return any(char.isdigit() for char in inputString)


##########################
## Main analysis window ##
##########################
#   I would like to add a user input: "blinded?"
#   I would like to add a 'save state' and 'load state' feature
#   to continue with previous analyses.     

i=0
while i < grand_total_traces: #switch to grand_total_traces for full version
    os.chdir(path)
    
    coordinates = get_coordinates(order[i])
    read_data(coordinates)
    plot()
    
    #options
    print(str(i+1)+"/"+str(grand_total_traces)) #progress marker
    ans = input('How many steps? Enter 1-8 or h for help: ')
    
    #forward, not enabled
    if ans == 'f':
        jump = input('How many do you want to skip? Enter an integer: ')
        i = i+int(jump)-1;
#
#    #go to, not enabled
    if ans == 'g': 
        print()
        for i in range(len(dir_names)):
            print(str(i)+" = "+dir_names[i])
        directory = int(input('Which directory? Enter a number: '))
        movie = int(input('Which movie? (starts from 1): '))
        trace = int(input('Which trace? (starts from 1): '))
        coordinates = [directory, movie, trace]
        read_data(coordinates)
        plot(False)
        i-=5; # I have no idea why this needs to be 5, else it skips ahead.
#
#    #save image not enabled
#    if ans=='s'
#        saveas(temp_image,'temp.jpg');  #JPG format
    
    if ans == 'h':
        print('''
              
              enter = pass / reject molecule
              b     = go back and reassign last
              x     = exit
              f     = skip forward
              g     = go to a molecule
              id    = display molecule coordinates
              
              In development:
              s     = save an image of the molecule
              
              ''')
        i-=1
    
    if ans == 'id':
        # I would like to find a way to not duplicate the plot when this happens.
        plot(False)
        i-=1
    
    #go back
    if ans == 'b':
        i-=1 # go back
        coordinates = get_coordinates(order[i]) # re-evaluate this file
        read_data(coordinates)
        os.chdir(current_analysis_path)
        os.chdir(dir_names[directory]+"/hel"+str(movie)) # go to the last folder
        tmp_parent = os.getcwd()
        find_name = "d"+str(directory)+"-m"+str(movie)+"-t"+str(trace) # trace to find
        #print('looking for '+find_name) # to match search terms
        step_dirs = sorted(next(os.walk(tmp_parent))[1]) # all the local directories
        for step_dir in step_dirs: # find last trace file and delete it
            os.chdir(step_dir)
            #print(step_dir) # to make sure you're looking everywhere
            for file in os.listdir():
                if file.split('_')[0] == find_name:
                    os.remove(file)
                    #print('Assignment for '+file+' deleted.') #unblinded only
                    print('Assignment for previous file deleted.')
                    print('Please reassign.')
            os.chdir(os.pardir)
        i-=1

    #enter, reject, pass, move on, etc. 
    elif ans in ('',0,'p'):
        i=i
        os.chdir(current_analysis_path)
        os.chdir(dir_names[directory]+"/hel"+str(movie)+"/Rejected")
        output_file_name = "d"+str(directory)+"-m"+str(movie)+"-t"+str(trace)+"_R.xlsx"
        data_list=[time, donor[trace-1], acceptor[trace-1]-leakage*donor[trace-1]]
        export(np.array(data_list),output_file_name) #saves a file
        os.chdir(path)
        
    # abort, exit, quit
    elif ans in ('x','q'):
        ok = input('Are you sure you want to quit? ')
        if ok in ('y',"ye","yes"):
            print('Analysis aborted by user.')
            break
            #sys.exit('Analysis aborted by user.') # if hard quit desired.
            

    #Photobleaching Steps
    #I would like to see if saving as a .csv or .txt would save space.
    elif has_numbers(ans):
        if int(ans) <= 8:
            step = int(ans)
            if step == 1:
                suffix = " step"
            else:
                suffix = " steps"
            os.chdir(current_analysis_path)
            os.chdir(dir_names[directory]+"/hel"+str(movie)+"/"+str(step)+suffix)
            output_file_name = "d"+str(directory)+"-m"+str(movie)+"-t"+str(trace)+"_"+str(step)+"s.xlsx"
            data_list=[time, donor[trace-1], acceptor[trace-1]-leakage*donor[trace-1]]
            export(np.array(data_list),output_file_name) #saves a file
            os.chdir(path)
            #j=j+1; #these counters were in the original matlab notebook.
            #mN(j)=i; #I haven't figured out what they're for.





    i+=1

##########################
#### Counting report #####
##########################
#   Here I would like to display a table with statistics
#   on the recently counted data. I would also like to 
#   make this a standalone script, so it must look nice!
#
#   I need to add an export function for this data.


os.chdir(current_analysis_path)
total_assigned = 0
total_rejected = 0
step_list = np.zeros([len(dir_names),8])

# Count
for directory in range(len(dir_names)):
    for movie in range(1,n_movies[directory]+1):
        os.chdir(dir_names[directory]+"/hel"+str(movie)+"/Rejected")
        total_rejected += len(os.listdir())
        os.chdir(current_analysis_path)
        for step in range(1,9):
            if step == 1:
                suffix = " step"
            else:
                suffix = " steps"
            os.chdir(dir_names[directory]+"/hel"+str(movie)+"/"+str(step)+suffix)
            total_assigned += len(os.listdir())
            step_list[directory,step-1] += len(os.listdir())
            os.chdir(current_analysis_path)

            
total=total_assigned+total_rejected
if total > 0:
    percent_assigned=(total_assigned/total)*100
else:
    percent_assigned=0
print()
print('total_assigned = '+str(total_assigned)+" ("+str(percent_assigned)+"%)")
print('total_rejected = '+str(total_rejected)+" ("+str(100-percent_assigned)+"%)")
print()

# print steps
header_list=['Steps']+[str(i)for i in list(range(1,9))]
print_list=[]
for i in range(len(dir_names)):
    print_list.append([dir_names[i]]+[j for j in step_list[i]])   
print(tabulate(print_list, headers=header_list))
print()

# print percentages
header_list=['Percentages']+[str(i)for i in list(range(1,9))]
print_list=[]
for i in range(len(dir_names)):
    print_list.append([dir_names[i]]+[(j/step_list[i].sum())*100 for j in step_list[i]])   
print(tabulate(print_list, headers=header_list))
    
