"""
File:           Organizer.py (Genome Sequencer Organization Script)

Author:         Alexander S. Adranly

Description:    Organizer.py is used to format and adjust the raw data that comes off the sequencer

"""

import os
import sys

###############################################
VERSION = 1.0
print "Genome Sequencer Organization: Version : " + str(VERSION)

SEQUENCER_DIRECTORY = '/Volumes/Hitachi 4/ACPHL Sequence Folder/031116'
COMPILED_DIRECTORY = 'COMPILED'

r2_store = []

###############################################

def _get_unique_sig(f):
    """
    :param: f: str : filename of the file to concatenate
    :return: str[ str (unique signature(without R#)), str[], str]

    description: creates a tuple with all the unique information about a specific file so that it can be downloaded
    """
    signature = f.split('_')
    return signature[0] + "_" + signature[1] + "_" + signature[2] + "_" + signature[4], signature, str(f)


def _is_sequenced_file(file):
    """
    :param file: str: the name of a file from the sequencer
    :return: boolean : if the file in question is a .fastq file

    description: Checks and confirms if the string in question is a sequenced filename
    """
    if file[-6:] == ".fastq":
        return True
    return False


def run():
    """
    :return:

    description: main function for this script
    """

    if not os.path.exists(os.path.expanduser(SEQUENCER_DIRECTORY)):
        print "Sequenced Directory Does not Exist!"
        return

    print "Valid Directory"

    if not os.path.exists(os.path.expanduser(SEQUENCER_DIRECTORY + '/' + COMPILED_DIRECTORY)):
        os.chdir(os.path.expanduser(SEQUENCER_DIRECTORY))
        os.mkdir(COMPILED_DIRECTORY)

    #set to sequence directory
    os.chdir(os.path.expanduser(SEQUENCER_DIRECTORY))
    directory_list = os.listdir(os.curdir)

    # setting up newly compiled files
    for raw_file in directory_list:

        if _is_sequenced_file(raw_file):
            # create a tuple of items based on its unique signature (void R#)
            unique_sig = _get_unique_sig(raw_file)
            signature = unique_sig[1]

            if signature[3] == "R1":
                # if R1, make a new file with its information in it
                try:
                    # open r1 to read from
                    in_file = open(raw_file)
                    # create a new file in the compiled directory with the unique sig
                    out_file = open(COMPILED_DIRECTORY + "/" + unique_sig[0], "w")
                    out_file.truncate()
                    # write contents of R1 into new file
                    out_file.write(in_file.read())
                except:
                    print "File I/O Error"
                finally:
                    in_file.close()
                    out_file.close()

            elif signature[3] == "R2":
                # if R2, store unique tuple into a tuple of r2 filenames
                r2_store.append(unique_sig)

    # concatenate the identical r2 to its newly compiled file
    try:
        os.chdir(COMPILED_DIRECTORY)
        comp_dir = os.listdir(os.curdir)

        for comp in comp_dir:
            for r2 in r2_store:
                if comp == r2[0]:
                    try:
                        in_file = open(SEQUENCER_DIRECTORY + "/" + r2[2])
                        out_file = open(comp, "a")
                        out_file.write("\n")
                        out_file.write(in_file.read())

                    except:
                        print "File I/O Error"
                    finally:
                        out_file.close()
                        in_file.close()

    except:
        print "Unmatched Compiler File Exception!" + str(sys.exc_info())
    finally:
        print "Genome Sequencer Organization Finished!"


if __name__ == '__main__':
    run()
