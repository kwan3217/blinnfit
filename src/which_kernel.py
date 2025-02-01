"""
Spice stuff, imported from emmexipy.python_lib.which_kernel
"""
from collections import namedtuple

import numpy as np
from spiceypy import kdata, cell_double, ckcov, spkcov, wncard, wnfetd, ktotal


ls_spice_return=namedtuple("ls_spice_return",("type","file"))
def ls_spice(verbose=False):
    """
    List all furnished spice kernels.
    :param verbose: if True, log also
    :return: A list of ls_spice_return named tuples.
    """
    count=ktotal('ALL')
    result=[]
    if verbose:
        print(f"Total of {count} kernels loaded")
    for i in range(count):
        (file,type,source,handle)=kdata(i, 'ALL')
        if verbose:
            print(f"{type:6s} {file}")
        result.append(ls_spice_return(type=type,file=file))
    return result


def which_kernel(type,obj,t,flush=False,print_cache=False):
    """
    Find out which loaded kernel provides a particular kind of data

    :param type:  one of 'LSK','FK','IK','SCLK' (collectively called text kernels)
                      or 'CK','SPK' (collectively called binary kernels)
    :param obj:   Ignored for text kernels, required for binary kernels,
                      numerical object id of spacecraft or frame to search for.
    :param t:     Ignored for text kernels, required for binary kernels,
                      Spice ET of time to search for. May be an array, see below.
    :param flush: if set, clear the cache.
    :return: path to kernel which covers this search, may be a list, see below.

    Notes

    - If looking for a text kernel, finds the last (highest priority) kernel of the given type, no matter what object it may have data for.
      Consequently, this only works well when data for only one spacecraft is loaded.

    - If looking for a binary kernel, you may pass an array of times. Return value is a list of kernels which provide data for this object
      at each time, or None if there is no kernel that covers the matching time. There is a one-to-one in-order relation between the list
      of times passed in and the list of files returned.

    - Return path is what was given in furnsh() (or in a furnished metakernel), and is either absolute or relative to what the current
      directory was at the time the kernel was furnished.

    - For performance reasons, this code keeps a cache in the form of function attribute. The cache is loaded whenever /flush is set, or
      the current cache is empty, or the requested object changes. Cache performs best therefore when you ask for coverage for the same object
      repeatedly, and switch objects a minimum number of times. You should pass /flush if you know that the kernels have changed since
      the last time you called which_kernel, because checking that the kernel list hasn't changed would take up too much time and invalidate
      a lot of the performance gain the cache provides.

    """
    utype=type.upper()
    if utype=='LSK' or utype=='FK' or utype=='IK' or utype=='SCLK':
        ftype='.T'+(utype[0:2] if len(utype)>2 else utype[0])
        count=ktotal('TEXT')
        for i in range(count-1,-1,-1):
            (this_file,this_type,source,handle)=kdata(i,'TEXT')
            #print,i,this_type,this_file,'   ',strlowcase(strmid(this_file,/rev,strlen(ftype)))
            if this_file[-len(ftype):].upper()==ftype:
                return this_file.replace('//','/')
    else:
        #Set up static variables
        if not hasattr(which_kernel,"windows"):
            which_kernel.windows={}
            which_kernel.files={}
        key=(utype,obj) #Python can index a dictionary with a tuple, so don't do things the IDL way (convert to string)
        if key in which_kernel.windows and not flush: #Check if the cache contains this object and data type
            #Retrieve the set of windows for this object and data type
            this_window=which_kernel.windows[key]
            this_files=which_kernel.files[utype]
        else:
            #Build the set of windows for this object and data type
            count=ktotal(utype)
            max_card=1
            #index 0 is file
            #index 1 is window in file
            #index 2 is begin (0) or end (1)
            #The effect is that this is a stack of matrices, where each plane is a file, each row in a plane is a
            #window, and each cell in the row is the beginning or end of the window. As the files are read, if a new
            #file has more windows than any previous file, exactly enough more rows are added to the bottom of all
            #the planes, and their value is set to NaN so they won't be found with np.where().
            this_window=np.zeros((count,max_card,2))*float('NaN')
            this_files=[""]*count
            for i_kernel in range(count):
                (this_file,this_type,source,handle)=kdata(i_kernel,utype)
                this_files[i_kernel]=this_file.replace('//','/')
                this_cover=cell_double(10000)
                if utype=='CK':
                    this_cover=ckcov(this_file, obj, False, 'INTERVAL', 0.0, 'TDB', this_cover)
                else:
                    this_cover=spkcov(this_file, obj, this_cover)
                card=wncard(this_cover)
                if card>max_card:
                    #Extend all the planes
                    this_window=np.concatenate((this_window,np.zeros((count,card-max_card,2))*float('NaN')),axis=1)
                for i_window in range(card):
                    (left,right)=wnfetd(this_cover,i_window)
                    this_window[i_kernel,i_window,0]=left
                    this_window[i_kernel,i_window,1]=right
            which_kernel.windows[key]=this_window
            which_kernel.files[utype]=this_files
        if print_cache:
            print(which_kernel.windows)
            print(which_kernel.files)
        try:
            n_elements_t=len(t)
        except:
            n_elements_t=1
        if n_elements_t==1:
            w = np.where(np.logical_and(t>=this_window[:, :, 0],t<=this_window[:, :, 1]))
            if len(w[0])>0:
                result=this_files[w[0][-1]]
            else:
                result=None
        else:
            result=[None]*n_elements_t
            for i_t in range(n_elements_t):
                w=np.where(np.logical_and(t[i_t]>=this_window[:,:,0],t[i_t]<=this_window[:,:,1]))
                if len(w[0])>0:
                    result[i_t]=this_files[w[0][-1]]
        return result
    raise RuntimeError('Somehow fell to the bottom of which_kernel')
