
# Ask indices by how much examples will be increased
try_again = True
while try_again:
    indices = input( "By how much do you want to increase the number of" +
                     "examples of each class ? \n For example: '1,2,10'" +
                     "does not change the number of mono_clustered, doubles" +
                     "the number of multi_clustered, and multiplies the " +
                     "number of EMFs by 10.")

    # Try to set the input string into a list of float
    try:
        indices = list( map( float, indices.split( ",")))

        # If ther is 3 items upper than 1 in the list, break the loop and continue
        if len( indices) == 3 and not any( x < 1 for x in indices):
            try_again = False
        # Otherwise, restart the loop
        else:
            print( "You entered too many indices, or one of them is below 1")
    # If the answer is another sting, restart the loop
    except:
        print( "Your answer is in a wrong format")


# Ask if the dataset generated will be in light_version
try_again = True
while try_again:
    light_version = input( "Light version of dataset ? (y/n): ")
    # If the answer is yes, stop the loop and create a variable suffix to load
    # light version of views
    if light_version == "y":
        suffix_light = "_light"
        try_again = False
    # If the answer is no, stop the loop and create a empty variable suffix to
    # load normal version of views
    elif light_version == "n":
        suffix_light = ""
        try_again = False
    # If the answer is another sting, restart the loop
    else:
        print( "You didn't enter \"y\" or \"n\"")

###############################################
# Import paths and modules
###############################################

# Define paths
DATA_DIR = "../../data/"
RAWDATA_DIR =  DATA_DIR + "rawData/"
H5_DIR = DATA_DIR + "dataset_compilation_to_hdf5/"

# Import modules
import pandas as pd
import numpy as np
import h5py
from imblearn.over_sampling import SMOTE, ADASYN

###############################################
# Oversampling
###############################################

# Load the dataset that will be oversampled
dataset = h5py.File( H5_DIR + "MuPPI_14view{}_EMF.hdf5".format( suffix_light), 'r')

# Create a list of view keys
views_keys = [ view for view in dataset.keys() if "View" in view]

# Count the number of proteins in each label class
labels_count = np.unique( dataset[ "Labels"][()], return_counts = True)[1]
# Multiply labels counts by indices entered
labels_count = np.array( labels_count * indices).astype( int)

# Store in a list the dataset view arrays
views = []
for key in views_keys:
    views.append( dataset[ key][()])

# concatenate those views by the columns
views = np.concatenate( views, axis = 1)
# Oversample the concatenated views with SMOTE
views_resampled, Labels_resampled = SMOTE({ 0:labels_count[ 0], 1:labels_count[ 1],
                                            2:labels_count[ 2]}).fit_resample(
                                            views, dataset[ "Labels"][()])


###############################################
# Set data in a new hdf5 file
###############################################

######## VIEWS ########

# Create a new hdf5 file for the oversampled data
dataset_resampled = h5py.File( H5_DIR + "MuPPI_14view{}_EMF_oversampled.hdf5", 'w')

# Define a start and stop integrer to select the features of the current
# view in the concatenated resampled views array. Set the start at 0
start = 0
# For each view key
for key in views_keys:
    # Set the stop equal to the number of columns of the current view plus the
    # start value
    stop = dataset[ key][()].shape[1] + start
    # Store the selected array in a hdf5 dataset
    dataset_resampled.create_dataset( key, ( views_resampled.shape[ 0], stop-start),
                                      data = views_resampled[:, start:stop])

    # In attribute of this dataset, stock the current view name, which is in
    # attribute in the initial dataset, and a boolean to say if the view is sparse
    dataset_resampled[ key].attrs[ "sparse"] = False
    dataset_resampled[ key].attrs[ "name"] = dataset[ key].attrs[ "name"]
    # Set the start variable at the current stop
    start = stop


######## LABELS ########

# Store in this hdf5 file the oversampled labels
dataset_resampled.create_dataset( "Labels", Labels_resampled.shape,
                                  data = Labels_resampled)
dataset_resampled[ "Labels"].attrs[ "names"] = dataset[ "Labels"].attrs[ "names"]


######## METADATA ########

# Create hdf5 file dataset "Metadata"
meta_data_grp = dataset_resampled.create_group( "Metadata")

# In attribute of this dataset, stock the number of view, number of class and
# number of exemples
meta_data_grp.attrs[ "nbView"] = len( views_keys)
meta_data_grp.attrs[ "nbClass"] = len( dataset[ "Labels"].attrs[ "names"])
meta_data_grp.attrs[ "datasetLength"] = len( Labels_resampled)

# Create a list of ids for created data
new_EMFs = []
for i in range( labels_count[ 1] - labels_count[ 2]):
        new_EMFs.append( "new_example_" + str( i))

# Convert the list in array
new_EMFs = np.array( new_EMFs).astype( np.dtype("S10"))
# Concatenate it with the proteins ids array of dataset
index_intesect = np.concatenate(( dataset[ "Metadata"][ "example_ids"][()],
                                  new_EMFs))

# In metadata groupe, create a dataset which contain the view indexes
meta_data_grp.create_dataset( "example_ids", index_intesect.shape,
                              data = index_intesect, dtype = np.dtype( "S10"))

# Close hdf5 files
dataset.close()
dataset_resampled.close()
