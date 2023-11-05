#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 11:35:16 2019

@author: Justin

----Counting Report----

Needs major clean-up

Generates a counting report from a path.
Must contain directories that contain movie direcotries (e.g. hel1, hel2, etc.)
Those directories must contain step subdirectories (1 step, 2 steps etc.) 
The files in those subdirectories are counted to make the final report.

Outline for this notebook:

1. Get the highest resolution data for steps / rejected in a master count list.
2. Format for movie and directory resolution exports
3. Print and export
"""

import os
#import glob
import sys
import numpy as np
from tabulate import tabulate
import pandas as pd

###############################
## Functions and definitions ##
###############################

# Path should be the path of an analysis performed, not the original data.
# Can add the functionality later to make it work on both
# Test path 1: /Users/Justin/Documents/Projects/InProgress/TM-EphA2-PIP2/Data/SMALPs/SM bleaching/2019-12-20 copy analyses/2019-12-23_07-51-05/
# Test path 2: /Users/Justin/Documents/Projects/InProgress/TM-EphA2-PIP2/Data/SMALPs/SM bleaching/2019-11-30/

path = "/Users/Justin/Documents/Projects/InProgress/TM-EphA2-PIP2/Data/SMALPs/SM bleaching/2020-02-20 analyses/2020-02-23_17-04-53/"
os.chdir(path);
dir_names = sorted(next(os.walk(path))[1]) # get directory names
n_dir = len(dir_names) # get the number of directories in path
print(n_dir," folders")


def get_n_movies(directory):
    '''
    input:
        int, directory index starting from 0
    output:
        returns the number of movies in the directory
    '''
    os.chdir(dir_names[directory])
    dirs = sorted(next(os.walk(os.getcwd()))[1])
    hel = 0
    for folder in dirs:
        if 'hel' in folder:
            hel+=1
    os.chdir(path)
    return hel


# Make a list of the number of movies in each directory.    
n_movies=[]    
for i in range(len(dir_names)):
    n_movies.append(get_n_movies(i))
print(n_movies," movies")


# get the nubmer of traces / molecules
def get_n_traces(directory, movie):
    '''
    input:
        int, directory index starting from 0
        int, movie index starting from 1
    output:
        returns the number of traces in a movie
    '''
    os.chdir(path)
    os.chdir(dir_names[directory])
    os.chdir("hel"+str(movie))
    subdirs = sorted(next(os.walk(os.getcwd()))[1])
    n_traces = 0
    for folder in subdirs:
        os.chdir(folder)
        counted = len([name for name in os.listdir('.') if os.path.isfile(name)])
        n_traces += counted
        os.chdir("..")
#    fname = "hel"+str(movie)+".traces"
#    fid = open(fname,'r') # the open file object  ### ERROR ###
#    length = int(np.fromfile(fid,'int32', 1)) # number of timesteps, must keep.
#    n_traces = int(np.fromfile(fid, 'int16', 1)) # number of traces
#    n_traces = int(n_traces/2)
#    fid.close()
    os.chdir(path)
    return n_traces


# New function
def get_step_dirs(directory, movie):
    '''
    Returns names of the '1 Step', '2 Steps' etc directories
    Only counts if contains 'step' or 'rejected'
    
    input:
        directory, int, index starting from 0
        movie, int, index starting from 1
    output:
        assign step_dirs with names of directories
    '''
    os.chdir(dir_names[directory]+"/hel"+str(movie))
    step_dirs = sorted(next(os.walk(os.getcwd()))[1])
    if len(step_dirs) == 0:
        sys.exit('''
                 Error: step directories not all present.
                 Data not fully analyzed.
                 Aborting.''')
    os.chdir(path)
    return step_dirs
step_dirs = get_step_dirs(0,1) # Assumes that all the directories are same.
  

# High resolution list for counting traces from each movie
os.chdir(path)
step_list = np.zeros([len(dir_names),len(step_dirs)])

# Gens empty list of architecture: 
# [slide1, [hel1, 0, 0, 0, 0, 0, 0, 0, 0],[hel2,..]...],[slide2...
master_count=[] 
for directory in range(len(dir_names)):
    master_count.append([dir_names[directory]])
    for movie in range(1,n_movies[directory]+1):
        master_count[directory].append(["hel"+str(movie)]+list(step_list[0]))

# Populate all lists with counts
for directory in range(len(dir_names)):
    for movie in range(1,n_movies[directory]+1):
        for step in range(len(step_dirs)):
            os.chdir(dir_names[directory]+"/hel"+str(movie)+"/"+step_dirs[step])
            master_count[directory][movie][step+1] += len(os.listdir())
            if 'Rejected' in step_dirs[step]:
                #total_rejected += len(os.listdir()) # Get rid of this
                step_list[directory,step] += len(os.listdir())
                os.chdir(path)
            else:
                #total_assigned += len(os.listdir()) # Get rid of this
                step_list[directory,step] += len(os.listdir())
                os.chdir(path)    


def calc_assigned(directory, movie):
    '''
    Function to calculate the total number of assigned molecules in a movie
    
    Input:
        directory, int, starting from zero
        movie, int, starting from one
    output:
        int, total number of molecules assigned in the movie
    '''
    sum_assigned = 0
    for i in master_count[directory][movie][1:-2]:
        sum_assigned += i
    return sum_assigned


def calc_rejected(directory, movie):
    '''
    Function to calculate the number of rejected molecules in a movie
    For use when rejected is not explicitly stated
    
    Input:
        directory, int, starting from zero
        movie, int, starting from one
    output:
        int, total number of molecules rejected in the movie
        Also, re-assigns rejected value to the number of rejected molecules.
    Dependencies:
        relies on master_count being defined
    '''
    global master_count
    if master_count[directory][movie][-1] == 0:
        sum_traces = 0
        for i in master_count[directory][movie][1:-1]:
            sum_traces += i
        n_traces = get_n_traces(directory,movie)
        n_rejected = n_traces - sum_traces
        master_count[directory][movie][-1] = n_rejected
    else:
        n_rejected = master_count[directory][movie][-1]
    return n_rejected


def recursive_sum(sum_function):
    '''
    Function to recursively apply a summation function to all subdirectories
    Input:
        sum_function, function, returns a sum from directory and movie
    Output:
        int, a sum of sums returned by the sum function
    '''
    tmp = 0
    for j in range(n_dir):
        for i in range(1,n_movies[j]+1): # check
            tmp += sum_function(j,i)
    return tmp


# Fix  rejected if needed
if step_list[0][-1] == 0:
    for directory in range(n_dir):
        for movie in range(1,n_movies[directory]+1):
            tmp = calc_rejected(directory, movie)
            step_list[directory][-1] += tmp

# Calculate totals
total_assigned = recursive_sum(calc_assigned)
total_rejected = recursive_sum(calc_rejected)

        
###############################
########## Printing ###########
###############################


def print_table(list_of_lists):
    '''
    Function to print large tables containing step counts for each movie
    
    Input:
        list_of_lists, list. First item is header.
        format: ['Slide 1 - 100-1',
                 ['hel1', 15.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 31.0],
                 ['hel2', 21.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 33.0]...]
    
    '''
    n_hel = len(list_of_lists)
    print_list=[]
    for i in range(1,n_hel):
        print_list.append([j for j in list_of_lists[i]]) 
    header_list=[list_of_lists[0]]+step_dirs
    print()
    print(tabulate(print_list, headers=header_list))


# Calculate and print totals
# I think this may be giving me total + 1 for some reason?            
total = total_assigned + total_rejected
if total > 0:
    percent_assigned = (total_assigned/total)*100
else:
    percent_assigned = 0

# Print statements of grand totals    
total_molecules_statement = 'total_molecules = '+str(total)
total_assigned_statement = 'total_assigned = '+str(total_assigned)+" ("+str(percent_assigned)+"%)"
total_rejected_statement = 'total_rejected = '+str(total_rejected)+" ("+str(100-percent_assigned)+"%)"  
print()
print(total_molecules_statement)
print(total_assigned_statement)
print(total_rejected_statement)
print()



# Print total steps and percentages for each Directory.
def print_totals_table(header,data_list):
    '''
    Prints a formatted table with headers
    Adds up totals and puts them at the end
    Assumes the last column is 'rejected'
    
    Input:
        header, str
        data, list of lists. no headers.
            Example: [[n1,n2...],[n3,n4...]...]
    Output:
        prints a table
    '''
    header_list = [header] + step_dirs + ['Total counted']
    print_list = []
    for i in range(len(dir_names)):
        tmp_count = [j for j in data_list[i]]
        sum_tmp_count = 0
        for k in tmp_count[0:-1]:
            sum_tmp_count += k
        print_list.append([dir_names[i]]+tmp_count+[sum_tmp_count])
    print(tabulate(print_list, headers=header_list))
    print()
    return [header_list,print_list]

# Ugly fix number one
percent_list = []
for i in range(len(step_list)):
    n_rejected = step_list[i][-1]
    n_total = step_list[i].sum()
    percent_rejected = (n_rejected/n_total)*100
    percent_list.append([(j/step_list[i][0:-1].sum())*100 for j in step_list[i][0:-1]]+[percent_rejected])
percent_list = np.array(percent_list)

total_steps_table = print_totals_table('Total steps',step_list)
total_percents_table = print_totals_table('Percentages',percent_list)

# Ugly fix number two
for i in range(len(total_percents_table[1])):
    total_percents_table[1][i][-1] = 100-total_percents_table[1][i][-2]
    

###############################
########### Export ############
###############################

def export(data_list,file_name="untitled.csv"):
    '''
    Exports data to file
    
    input:
        data_list, list of lists of ints or floats. Must be ndarray.
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
    #df = df.T # transpose the data
    df.to_csv(file_name,index=None,header=None)

# Export individual counts
def format_for_export(list_of_lists):
    header_list=[list_of_lists[0]]+step_dirs
    array = np.asarray(list_of_lists[1:],dtype='<U50')
    array = np.insert(array,0,header_list,axis=0)
    array = np.insert(array,len(list_of_lists[1:])+1,'',axis=0)
    return array
    
data_list = []
for i in range(len(master_count)):    
    tmp_list = format_for_export(master_count[i])
    data_list.append(tmp_list)
data_list = np.concatenate(data_list)

directory_basename = os.path.basename(os.path.normpath(path))
out_file_name = directory_basename+' individual_counts.csv'
export(data_list,path+out_file_name)

# Export totals and percentages
# Appending one row at a time
# Not an elegant solution. I know.

n_cols = len(total_steps_table[0])

export_list = np.asarray(total_steps_table[1],dtype='<U50')
export_list = np.insert(export_list,0,total_steps_table[0],axis=0) # first headers
export_list = np.insert(export_list,len(export_list),'',axis=0) #from here on adding rows to the end one by one.
export_list = np.insert(export_list,len(export_list),total_percents_table[0],axis=0)
export_list = np.insert(export_list,len(export_list),total_percents_table[1],axis=0)
export_list = np.insert(export_list,len(export_list),'',axis=0)
export_list = np.insert(export_list,len(export_list),[total_molecules_statement]+['']*(n_cols-1),axis=0)
export_list = np.insert(export_list,len(export_list),[total_assigned_statement]+['']*(n_cols-1),axis=0)
export_list = np.insert(export_list,len(export_list),[total_rejected_statement]+['']*(n_cols-1),axis=0)
export_list = np.insert(export_list,len(export_list),'',axis=0)
export_list = np.insert(export_list,len(export_list),['Generated by Sum_counted.py']+['']*(n_cols-1),axis=0)


out_file_name = directory_basename+' count_totals.csv'
export(export_list,path+out_file_name)



