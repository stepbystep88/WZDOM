

Bio_ and Phy_ Datasets Format:   http://kodiak.cs.cornell.edu/kddcup/datasets.html
===================================================================================


    * phy_train.dat: Training data for the quantum physics task (50,000 train cases)
    * phy_test.dat: Test data for the quantum physics task (100,000 test cases)
    * bio_train.dat: Training data for the protein homology task (145,751 lines)
    * bio_test.dat: Test data for the protein homology task (139,658 lines) 

The file formats for the two tasks are as follows.
Format of the Quantum Physics Dataset
Each line in the training and the test file describes one example. The structure of each line is as follows:

    * The first element of each line is an EXAMPLE ID that uniquely describes the example. You will need this EXAMPLE ID when you submit results.
    * The second element is the class of the example. Positive examples are denoted by 1, negative examples by 0. Test examples have a "?" in this position. This is a balanced problem so the target values are roughly half 0's and 1's.
    * All following elements are feature values. There are 78 feature values in each line.
    * Missing values: columns 22,23,24 and 46,47,48 use a value of "999" to denote "not available", and columns 31 and 57 use "9999" to denote "not available". These are the column numbers in the data tables starting with 1 for the first column (the case ID numbers). If you remove the first two columns (the case ID numbers and the targets), and start numbering the columns at the first attribute, these are attributes 20,21,22, and 44,45,46, and 29 and 55, respectively. You may treat missing values any way you want, including coding them as a unique value, imputing missing values, using learning methods that can handle missing values, ignoring these attributes, etc. 

The elements in each line are separated by whitespace.
Format of the Protein Homology Dataset
Each line in the training and the test file describes one example. The structure of each line is as follows:

    * The first element of each line is a BLOCK ID that denotes to which native sequence this example belongs. There is a unique BLOCK ID for each native sequence. BLOCK IDs are integers running from 1 to 303 (one for each native sequence, i.e. for each query). BLOCK IDs were assigned before the blocks were split into the train and test sets, so they do not run consecutively in either file.
    * The second element of each line is an EXAMPLE ID that uniquely describes the example. You will need this EXAMPLE ID and the BLOCK ID when you submit results.
    * The third element is the class of the example. Proteins that are homologous to the native sequence are denoted by 1, non-homologous proteins (i.e. decoys) by 0. Test examples have a "?" in this position.
    * All following elements are feature values. There are 74 feature values in each line. The features describe the match (e.g. the score of a sequence alignment) between the native protein sequence and the sequence that is tested for homology.
    * There are no missing values (that we know of) in the protein data. 

To give an example, the first line in bio_train.dat looks like this:

279 261532 0 52.00 32.69 ... -0.350 0.26 0.76

279 is the BLOCK ID. 261532 is the EXAMPLE ID. The "0" in the third column is the target value. This indicates that this protein is not homologous to the native sequence (it is a decoy). If this protein was homologous the target would be "1". Columns 4-77 are the input attributes.
The elements in each line are separated by whitespace.