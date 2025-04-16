"""This file defines useful functions for dataset generation scripts."""

# Import paths
RAWDATA_DIR = "../../data/rawData/"

# Import modules
import pandas as pd
import numpy as np


#########################################
# Functions
#########################################

def explodeValues(df, index_col, col_to_explode, separator):
    """This function separates the str values of a column into several rows. It
    takes as input a df with a column as index and a column whose values are to
    be separated"""

    # Replace NaN by empty strings
    df = df.replace(np.nan, '', regex=True)
    # Create a new df with as index index_col and as single column
    # col_to_explode whose values are separated by the separator
    df = pd.DataFrame(df[col_to_explode].str.split(separator).tolist(),
                      index=df[index_col]).stack()
    # Turn the indexes into a column again
    df = df.reset_index([0, index_col])
    df.columns = [index_col, col_to_explode]

    # Return this two columns df
    return (df)


def createTable(df, col_to_indexes, col_to_features):
    """This function creates a df of 0 and 1 from two df columns, one converted
    to an index, the other to column names. The col_to parameters must be str"""

    # Keep only the col_to_indexes and col_to_features columns in df
    df = df[[col_to_indexes, col_to_features]]
    # Create a column withs stinrgs that concatenate the two df columns
    df = df.assign(pair=df[col_to_indexes] + df[col_to_features])
    # Count the number of those strings and put them in a new column
    val = df.pair.value_counts().reset_index()
    df = df.drop_duplicates()

    # CHANGED: Using the right col name and correct the merge
    val['count'] = val['count'].astype(str)
    df = df.merge(val, on="pair")
    # Pivot the df with col_to_indexes as index, col_to_features as columns,
    # and the value counts as values
    df = df.pivot(index=col_to_indexes, columns=col_to_features,
                  values='count')
    # Convert np.nan to zeroes
    df = df.replace(np.nan, 0).astype(int)

    return (df)


def fillAncestorsTerms(table, heir_term, anc_list):
    """This function adds to an ontology table the features of an heir term's
    ancestor terms, and adds them together"""

    # For ancestors already present: we add the heir_term column
    present_anc = [anc for anc in anc_list if anc in table.columns]
    for anc in present_anc:
        table[anc] += table[heir_term]

    # For ancestors that do not yet exist: we prepare a DataFrame
    new_anc = [anc for anc in anc_list if anc not in table.columns]
    if new_anc:
        # Create a new dataframe with the columns
        new_cols = pd.DataFrame({anc: table[heir_term] for anc in new_anc}, index=table.index)
        # CHANGED: Use of pd.concat
        table = pd.concat([table, new_cols], axis=1)

    return table


def uniquePairs(arr):
    """ This function identifies the identical lines of a numpy array, and
    returns a list of lists with in each one the indexes of the identical lines"""

    # Converts the values of each row into a single variable with all values inside.
    uview = np.ascontiguousarray(arr).view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))
    # Stores the unique values in a array uvals, and their indexes in the matrix
    # in a array uidx
    uvals, uidx = np.unique(uview, return_inverse=True)
    # Searches uidx for uvals indexes that are present more than once, and
    # stores them in another array
    pos = np.where(np.bincount(uidx) > 1)[0]
    # For each index duplicated in uidx (and therefore identical lines in the
    # matrix), look for the position of these indexes in uidx, which corresponds
    # to the position of the identical lines in the initial matrix. Stores the
    # list of these positions in an array which is returned as output.
    pairs = []
    for p in pos:
        indices = np.where(uidx == p)[0]
        pairs.append(indices.tolist())

    return pairs


def delRedundantCols(table):
    """ This function identifies the identical columns of a pandas DataFrame,
    removes redundant columns and names the remaining one as a concatenation of
    all identical column names"""

    # Stores the names of identical columns, i.e. terms that are annotated for
    # exactly the same proteins
    duplicated_cols = table.columns[table.T.duplicated(keep=False)]
    identical_col_lists = uniquePairs(table[duplicated_cols].values.T)

    # Iterate the found pairs, store in a key dict the name of the column
    # corresponding to the first index of the pair, then in value the names of
    # the concatenated columns of all the indexes of the list of pairs.
    identical_col_dict = {}

    for identical_col in identical_col_lists:
        identical_col_dict[duplicated_cols[identical_col[0]]] = '_'.join(
            duplicated_cols[identical_col])

    # Lists duplicate collars that are not in the keys of the dict and removes
    # them from the df
    col_to_delete = list(set(duplicated_cols) - set(identical_col_dict.keys()))
    table = table.drop(col_to_delete, axis=1)

    # Renames the remaining duplicate cols (those in dict keys) by the str of
    # all identical cols (the dict values).
    table = table.rename(identical_col_dict, axis=1)
    return (table)


def CutSmallLeafs(df, leaf_max_size=6):
    """This feature reduces the number of features in an ontology Boolean
    DataFrame. It removes features that annotate less than N proteins. The
    number N is set by the parameter \"leaf_max_size\", which defaults to 5."""

    # Convert df to np array of 0 and 1
    arr = df.values
    arr = np.where(arr > 0, 1, 0)
    # Identifies columns whose sum of 1's is greater than or equal to the
    # leaf_max_size parameter, and returns the df with only those columns.
    mask = np.where(sum(arr) >= leaf_max_size, True, False)
    return (df[df.columns[mask]].astype(int))


def MissingExamples(df, zero_or_nan):
    """This function take a pandas DataFrame as entry, with protein ids as index.
    It add the missing PPI network protein ids in the index, and fill the columns
    with zeroes or numpy nan"""

    # Load the list of PPI proteins
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
        PPInetwork_proteins = f.read().split('\n')
        PPInetwork_proteins.remove('')

    # List the missing proteins in the df indexes
    missed_prot = [prot for prot in PPInetwork_proteins if
                   prot not in df.index]

    # Create a new df with the columns of the input df and the missing
    # proteins as an index. If the zero_or_nan argument is a "zero" string
    # fill the rows with zeroes
    if zero_or_nan == "zero":
        df_MissEx = pd.DataFrame(np.zeros((len(missed_prot), len(df.columns))),
                                 columns=df.columns, index=missed_prot)
    #  If the zero_or_nan argument is a "nan" string,fill the rows with numpy nan
    elif zero_or_nan == "nan":
        arr = np.empty((len(missed_prot), len(df.columns)))
        arr[:] = np.nan
        df_MissEx = pd.DataFrame(arr, columns=df.columns, index=missed_prot)

    # Concatenate the two df and return it
    df = pd.concat([df, df_MissEx])
    return (df)
