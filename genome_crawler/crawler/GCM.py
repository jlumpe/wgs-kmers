"""
File:           GCM.py (Genome Compilation Module)

Author:         Alexander S. Adranly

Description:    GCM.py contains classes and functions that are elements of the system necessary to
                all of the requested bacteria genome files from ftp.ncbi.nlm.nih.gov

classes:        GCM:            Object in charge of creating, orchestrating, and allocating resources to the spider and
                                filemanager
                Spider:         Downloads the files targeted towards it
                TargetManager:  Allocates and processes the filepaths for the spider to more easily download everything
                Filemanager:    File preparation for the spider to download the needed files


NOTE:
    deal with errors in the individual thread instead of working on files back and forth
    have the thread restart itself
"""
import time
import ftputil
import threading
import os
import sys

#######################################################
# Global Variables
VERSION = 1.5
DEBUG = True
print "Genome Compilation Module: Version: ", VERSION, '\n'

# DIRECTORIES ON FTP SERVER
FTP_HOST_PATH = 'ftp.ncbi.nlm.nih.gov'  # FTP HOST
BACTERIA_PATH = '/genomes/refseq/bacteria/'
ALL_GENOME_PATH = '/genomes/all/'
LATEST_PATH = '/latest_assembly_versions/'
ALL_PATH = '/all_assembly_versions/'

# Filepath for destination of downloaded files
TARGET_DIRECTORY = "GENOME_DATA"  # "/Users/Silversmith/Desktop/DATA/"
THUNDER_DIRECTORY = "/Volumes/Hitachi 4/"
GENERIC_DIRECTORY = "~/Desktop/"

CHOICE_DIRECTORY = GENERIC_DIRECTORY

# Filepath integrity Tuples
NON_DIVEABLES = ".txt"
DIVEABLES = "GCF"
DOWNLOADABLES = ".fna.gz"

# Threading Controls
THREAD_COUNT = 10
SLEEP_TIME = 30  # time also to wait after a connection has been rejected
#######################################################

"""
class FileManager(object)

:An object in charge of organizing/creating/and destroying files (file preparation) before downloads start
NOTE:
Will be adding a feature that finds out what the genus and species is... hopefully
this may also possibly go in the spider section since the spider really goes deep in directory
"""


class FileManager(object):
    def _create_base_directory(self):

        if not os.path.exists(os.path.expanduser(CHOICE_DIRECTORY + TARGET_DIRECTORY)):
            os.chdir(os.path.expanduser(CHOICE_DIRECTORY))
            os.mkdir(TARGET_DIRECTORY)
        else:
            print "Directory Exists"

    def bacteria_setup(self, bacteria_list):
        """
        :param bacteria_list: str[] : list of bacteria directories
        :return: void

        description: creates a directory group for storing the downloaded files
        """
        print "\n ___________________________________________"
        print "Creating Download Directory (~/Desktop/GENOME_DATA)..."
        self._create_base_directory()
        os.chdir(os.path.expanduser(CHOICE_DIRECTORY + TARGET_DIRECTORY))
        print "Creating Relevent Genome Directories..."

        if len(bacteria_list) == 0:
            print "Directory Setup Complete!"
            print " ___________________________________________\n"
            return

        for bacteria in bacteria_list:
            if not os.path.exists(os.path.expanduser(CHOICE_DIRECTORY + TARGET_DIRECTORY + '/' + bacteria)):
                os.mkdir(bacteria)

        print "Directory Setup Complete!"
        print " ___________________________________________\n"


"""
class Spider(threading.Thread):

:A class responsible for downloading all of the bacteria file targets distributed to it from the TargetManagers
"""


class Spider(threading.Thread):
    def __init__(self, gcm):
        threading.Thread.__init__(self)
        self.GCM = gcm
        self.target_list = []
        self.host = None

        # flags
        self.target_files_downloaded = 0

    def get_remaining_targets(self):
        """
        :return:
        """
        return self.target_list

    def _is_downloadable(self, file_name):
        """
        :param file_name: str :name of file or directory

        :return: boolean: returns if this is a downloadable .fna file

        description: confirms or file downloadability given a file
        """

        if len(file_name) > 24:
            if file_name[-24:] == "_cds_from_genomic.fna.gz":
                return False

        if len(file_name) > 24:
            if file_name[-24:] == "_rna_from_genomic.fna.gz":
                return False

        if len(file_name) > 15:
            if file_name[-15:] == "_genomic.fna.gz":
                return True

        return False

    def _download(self, file_list, name):
        """
        :param file_list: str[] : list of str elements, each the name of a file in the directory
        :param name: str : name of bacteria downloaded

        descrition: given a directory and a file name, the appropriate files are downloaded
        """
        for file in file_list:
            try:
                if self._is_downloadable(file):
                    if not DEBUG:
                        os.chdir(os.path.expanduser(CHOICE_DIRECTORY))
                        self.host.download(file, TARGET_DIRECTORY + '/' + name + '/' + file)
                    print "@" + self.getName() + " " + name + " : " + file  # print statement
                    self.target_files_downloaded += 1

            except OSError:
                print "Directory Does Not/No Longer Exists!"


    def _get_name(self, path):
        """
        :param path: str: GCF directory path that holds .fna files
        :return: str : name of the bacteria of this GCF path

        description: given a filepath, the bacteria name is obtained
        """
        path_parts = path.split('/')
        return path_parts[4]

    def _connect_host(self):
        """
        :return:

        Description: Method used to connect host to server
        """
        self.host = None  # initialize it to None, so we dont deal with false positives

        while self.host is None:
            print "@" + self.getName() + ": Attempting to Connect Host..."
            try:
                self.host = ftputil.FTPHost(FTP_HOST_PATH, 'anonymous', 'password')  # personal host ftp connection
            except OSError:
                print "@" + self.getName() + ": Error While Attempting to Connect Host!"
                time.sleep(SLEEP_TIME)

        print "@" + self.getName() + ": Host Connected!"

    def crawl(self):
        """
        :return: void

        description: goes through all targeted GCF files and downloads all of them
        """

        while len(self.target_list) > 0:
            path = self.target_list.pop()

            try:
                # Assuming that you are given a GCF path
                self.host.chdir(path + "/")  # enter path
                self._download(self.host.listdir(self.host.curdir), self._get_name(path))  # download .fna file

            except OSError:
                print "@" + self.getName() + ": OS Error thrown and Caught" + str(sys.exc_info())
                print "\tConnection Lost..."
                self.host.close()  # possibly this may stop the timeout errors from occuring
                self.target_list.append(path)  # re-add the element that the error occured in
                time.sleep(SLEEP_TIME)
                self._connect_host()

            except:
                print "@" + self.getName() + "General Exception" + str(sys.exc_info())

    def run(self):
        """
        :return: void

        Description: starts when the start() method is called
        """

        if len(self.target_list) == 0:
            print self.getName() + ": No targets assigned to Spider, CLOSING!"
            return

        self._connect_host()
        self.crawl()
        self.host.close()  # must finally close host

"""
class TargetManager(threading.Thread):

:A class responsible for extracting all the targets for the spider to download
"""


class TargetManager(threading.Thread):
    def __init__(self, gcm):
        threading.Thread.__init__(self)
        self.GCM = gcm
        self.directory_list = []  # stores the bacteria files to go through (init)
        self._expanded_list = []  # stores the processed directories for the spiders to crawl (fin)
        self.host = None

    def get_expansion(self):
        """
        :return: str[]

        Description: returns the list of GCF directories to download fna files from
        """
        return self._expanded_list

    def is_suppressed(self, item):
        """
        :param item: str: the folder which the code intends to go down
        :return: bool: if this folder is labeled 'suppressed'

        Description: Returns a boolean verifying if the string matches a specific phrase
        """
        return item == "suppressed"

    def _has_assembly(self, assembly_list):
        """
        :param assembly_list:
        :return: tuple(Boolean, Boolean) : stores if there is a latest or all assembly directory

        Description: Checks if the directory passed has the assembly folders needed to continue
        """
        latest = False
        all = False
        for assembly in assembly_list:
            if assembly == "latest_assembly_versions":
                latest = True
            if assembly == "all_assembly_versions":
                all = True

        return latest, all

    def _is_diveable(self, dir_name):
        """
        :param dir_name: str : the directory in question
        :return: boolean: if the file is a direcotry we should go into

        Description: Method used to decide if a file/directory is "diveable"
                     (do we want to go inside the directory?)
        """
        if dir_name[-4:] == NON_DIVEABLES:
            return False

        return True

    def _connect_host(self):
        """
        :return:

        Description: Method used to connect host to server
                     If it fails the first time, it will continue to try
        """
        self.host = None  # initialize it to None, so we dont deal with false positives

        while self.host is None:
            print "@" + self.getName() + ": Attempting to Connect Host..."
            try:
                self.host = ftputil.FTPHost(FTP_HOST_PATH, 'anonymous', 'password')  # personal host ftp connection
            except OSError:
                print "@" + self.getName() + ": Error While Attempting to Connect Host!----\n"
                time.sleep(SLEEP_TIME)

        print "@" + self.getName() + ": Host Connected!"

    def _expand_directory_list(self):
        """
        :param list:
        :return: list : list of all GCFs needed to downlad

        description: method that takes the directory of bacteria and makes a list off all the GCFs to look into
        """
        latest = False  # stores if the current directory has the latest assembly
        expanded_list = []  # for storing all files

        while len(self.directory_list) > 0:

            element = self.directory_list.pop()
            # for element in self.directory_list:
            #print element  # supressing print

            # initial filtering
            if not self._is_diveable(element):
                continue

            # 1. Go into Directory name
            modified_path = BACTERIA_PATH + element

            try:
                self.host.chdir(modified_path)  # Entering the Directory

                # look for latest assembly and all assembly
                _assemblies = self._has_assembly(self.host.listdir(self.host.curdir))

                if _assemblies[0]:
                    # if has Latest
                    modified_path += LATEST_PATH  # adding latest assembly onto path
                    latest = True
                elif _assemblies[1]:
                    # if has all
                    modified_path += ALL_PATH  # adding latest assembly onto path
                    latest = False
                else:
                    # if has none
                    continue

                # move into assembly
                self.host.chdir(modified_path)
                gcf_list = self.host.listdir(self.host.curdir)  # get list of GCFs
                # there may be suppressed folders inside!!

                for gcf_file in gcf_list:
                    # iterate through each gcf file

                    if self.is_suppressed(gcf_file):
                        # print element + ": has a suppressed folder!"
                        continue  # skip all the files in here

                    modified_path += gcf_file
                    expanded_list.append(modified_path)

                    # reset modified path after every iteration so there aren't any snowball effects
                    if latest:
                        modified_path = BACTERIA_PATH + element + LATEST_PATH
                    else:
                        modified_path = BACTERIA_PATH + element + ALL_PATH

            except OSError:
                print "@" + self.getName() + ": OS Error thrown and Caught" + str(sys.exc_info())
                print "\tConnection Lost..."
                self.host.close()  # possibly this may stop the timeout errors from occuring
                self.directory_list.append(element)  # re-add the element that the error occured in
                time.sleep(SLEEP_TIME)
                self._connect_host()

            except:
                print "@" + self.getName() + "General Exception" + str(sys.exc_info())

        # save personally so list can be used later
        self._expanded_list = expanded_list

    def run(self):
        """
        :return: void

        Description: starts when the start() method is called
        """
        if len(self.directory_list) == 0:
            print self.name + ": No targets assigned to Manager!"
            print "Closing: " + self.name
            return

        self._connect_host()
        self._expand_directory_list()  # allocate targets for spiders
        self.host.close()  # must finally close host



"""
class GCM(object):

:The object in charge of scheduling and allocating targets for the spiders to download
Main host for the orchestration of file downloading and file organization
"""


class GCM(object):
    def __init__(self):
        self.host = None
        # for assurance of Accuracy
        self.target_files_to_download = 0
        self.target_files_downloaded = 0

        self._spider_list = []  # stores spiders
        self._manager_list = []  # stores managers
        self.TARGET_LIST = []

        # create tuple of managers and spiders
        for i in xrange(0, THREAD_COUNT):
            self._manager_list.append(TargetManager(self))
            self._spider_list.append(Spider(self))

    def _generate_time_estimate(self):
        """
        :return:
        """
        print "I hope to estimate how long it takes sometime"

    # reporter function for results
    def _generate_report(self):
        print "\tREPORT:\n"
        print "Number of Threads: ", THREAD_COUNT
        print "Number of Target Files To Download: ", self.target_files_to_download
        print "Number of Target Files Downloaded: ", self.target_files_downloaded
        for spider in self._spider_list:
            print spider.getName() + " Download Count: " + str(spider.target_files_downloaded)

    def _clean(self):
        self.target_files_to_download = 0
        self.target_files_downloaded = 0

        # reset the spider list
        self._spider_list = []
        self._manager_list = []
        for i in xrange(0, THREAD_COUNT):
            self._manager_list.append(TargetManager(self))
            self._spider_list.append(Spider(self))

        self.thread_error_queue = []
        self.TARGET_LIST = []

    #  crawling internals
    def _generate_target(self, dir_list):
        """
        :param dir_list: str[] : a list of filepath information to distribute
        :return: item : an indexed item from an iterable object

        description: a generator to distribute any information necessary
        """
        for dir in dir_list:
            yield dir

    def _target_manager_sequence(self, dir_list):
        """
        :param dir_list:
        :return:

         description: A method that finds all the target GCF files
        """

        # TargetManager Routine

        dir_gen = self._generate_target(dir_list)  # generator for managers

        try:
            while True:
                for manager in self._manager_list:
                    manager.directory_list.append(next(dir_gen))

        except StopIteration:
            print "\n***Target Distribution Complete***\n"

        # run all managers
        for manager in self._manager_list:
            manager.start()
            time.sleep(1)  # just to let the manager before have a head start

        for manager in self._manager_list:
            manager.join()

        for manager in self._manager_list:
            self.TARGET_LIST += manager.get_expansion()  # whats already been parsed

        self.target_files_to_download += len(self.TARGET_LIST)  # it just keeps on adding up

    def _spider_sequence(self):
        """
        :return:
        """
        # SPIDERS Routine
        # generator for Spiders

        target_gen = self._generate_target(self.TARGET_LIST)

        try:
            while True:
                for spider in self._spider_list:
                    spider.target_list.append(next(target_gen))

        except StopIteration:
            print "\n***Spider Distribution Complete***\n"

        # run all spiders
        for spider in self._spider_list:
            spider.start()

        for spider in self._spider_list:
            spider.join()

        # Collect all download information for report generation
        for spider in self._spider_list:
            self.target_files_downloaded += spider.target_files_downloaded

    def crawl_genomes(self, bacteria_filename=None):
        """
        :param bacteria_filename: str|None
        :return: None

        description: method that finds all the directories to download and allocates them to spiders
        """
        print "\n ___________________________________________\n CRAWLING BACTERIA"

        manager = FileManager()  # for file organization

        try:
            self.host = ftputil.FTPHost(FTP_HOST_PATH, 'anonymous', 'password')  # attempt connection with ftp
            directory_list = []  # initialize the object that stores all the targets

            # for download query requests
            if bacteria_filename is None:
                # if there is not a query request
                self.host.chdir(BACTERIA_PATH)  # Changing directory to ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria
                directory_list = self.host.listdir(
                    self.host.curdir)  # Retrieving list of directories and storing it in a list
                print str(len(directory_list)) + " Directories"

            elif bacteria_filename is not None:
                # if there is a query request
                directory_list = bacteria_filename
                print str(len(directory_list)) + " Directories"
                print directory_list

            manager.bacteria_setup(directory_list)  # has to be done before any downloading takes place
            self.host.close()  # close ftp connection
        except OSError:
            print "Connection Attempt Failed"
            return

        try:
            # subdivide all of the bacteria files into their GCF files
            print "\n ****************************\n\tMANAGER SEQUENCE\n****************************\n"
            self._target_manager_sequence(directory_list)  # targetmanager
            print "\n ****************************\n\tSPIDER SEQUENCE\n****************************\n"
            self._spider_sequence()
        finally:
            print "\nCrawling Finished!\n"
            print " ___________________________________________\n"
            self._generate_report()
            print " ___________________________________________\n"
            self._clean()
