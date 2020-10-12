# This program is a user interactive program that parses PDB files that the user provides it with. 

# All modules, packages and libraries that are used in the program are imported.

import os.path
import sys
import os
import time
import numpy
from Bio.PDB.PDBParser import PDBParser
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from matplotlib import style


# This function will be called when the user input limit is reached in function_1(), function_2() and function_4().

def Redirection_to_change_directory():

    working_directory_path = os.getcwd()

    print("\n\nYour current working directory is: " + working_directory_path)

    directory_path()
    Main_Menu()   


# The function function_1() takes a PDB file from the user and outputs structural infromation.
# The structural infromation outputted to the user is the number of ATOM and HETATM entries of the protein and how many unique chains the protein has.
# The variable input_limit_counter_function_1 is initialised to 0.
# This variable is there to make sure that once the user enters a consecutive number of invalid inputs, they are redirected to allow them to change their working directory.
# This makes sure that the user knows what their current working directory for the program is and to prompt them to change it.
# This is to make sure that the directory that the user selects contains PDB files that can be used in the program.
# The while loop contains if conditional statements checking whether the user has inputted a PDB file and if the PDB os.path.exists.
# If the user enters a PDB file that exists in the working directory then that file will be opened.
# The ATOM_count and HETATM_count variables are intialised to 0.
# The for loop iterates over all the lines of the PDB file provided by the user.
# The startswith() method is used so that if a line in the PDB file starts with the string "ATOM", the ATOM_count increases by one.
# The startswith() method is also used so that if a line in the PDB file starts with the string "HETATM", the HETATM_count increases by one.
# The third if statement in the for loop uses the startswith() method to find lines that start with the string "COMPND".
# In the same if statement, the and is used so the line has to both start with "COMPND" and contain the string "CHAIN:".
# The replace_string variable replaces the string and prints the new string along with parts of the line that were not replaced which are the unique chains of the protein.
# The user is then prompted and asked for input on whether they want to return to the Main Menu or exit the program.
# The while loop is contructed so that the users input can be validated.
# If the user enters MENU, the Main_Menu() function is called and the program returns to the Main Menu.
# If the user enters EXIT, the exit() method stops the program running.

def function_1():

    input_limit_counter_function_1 = 0

    user_entry = input("Please enter the PDB file you would like to work with: ")

    while True:

        if input_limit_counter_function_1 == 5:
            print("\nYou have not entered a PDB file within your working directory\n")
            print(
            "You will be redirected to the main page\n"
            + "This will allow you to change your directory as you probably have no PDB files in this directory\n\n")
            time.sleep(3)
            Redirection_to_change_directory()
            break
            
        if user_entry [-3:] != "pdb":
            input_limit_counter_function_1 +=1
            user_entry = input("This is not a PDB file, please enter a pdb file: ")
            continue

        if not os.path.exists(user_entry):
            input_limit_counter_function_1 += 1
            print("The PDB file %s does not exist in your directory"%user_entry)
            user_entry = input('Please try again: ')
            continue
    
        if user_entry[-3:] == "pdb" and os.path.exists:

            PDB_File = open(user_entry)

            ATOM_count = 0
            HETATM_count = 0
    


            for line in PDB_File:
                if line.startswith("ATOM"):
                    ATOM_count += 1 
                if line.startswith("HETATM"):
                    HETATM_count += 1
                if line.startswith("COMPND") and line.__contains__("CHAIN:"):
                    replace_string = line.replace("COMPND   3 CHAIN:","This unique chains in this protein are: ")
                    print("\n" + replace_string)
         
        
            print("Number of ATOM entries in the protein: " + str(ATOM_count))
            print("Number of HETATM entires in the protein: " + str(HETATM_count))
            PDB_File.close()
            break

            
            
    menu_or_exit = input("\nEnter MENU to return back to the Main Menu and EXIT to exit the program\n")
    
    while True: 
        if menu_or_exit == "MENU":
            Main_Menu()
            break

        if menu_or_exit == "EXIT":
            exit()

        if menu_or_exit != "MENU" or menu_or_exit != "EXIT":
            menu_or_exit = input("\nThat is not a valid option!" + 
            "\nEnter Menu to return back to the menu page and EXIT to exit the program\n")


# The function function_2() calculates C-alpha distance between two C-alpha atoms of two different residues in the protein.
# The variable input_limit_counter_function_2 is intialised to 0 and this will make sure there is a limit on the number of invalid inputs from the user.
# This is done to remind the user that they might not have any PDB files in the current working directory of the program.
# This is so they can be redirected to change the working directory to one where PDB files are present that they can work with in the program.
# The input prompt asks the user to enter the filename of the PDB file they want to work with.
# The while loop has if conditional statements and a else statement, these check whether the user has entered a PDB file and if os.path.exists or not.
# If the filename that the user has entered does not end with the .pdb extension then the user will be told that their input is not a pdb file.
# If os.path does not exist then the user will be told that the PDB file they have entered does not exist in the working directory they are using for the program.
# If the user has entered a PDB file in the working directory and the os.path.exists then C-alpha distance can be calculated.
# Using the PDBParser from Bio.PDB.PDBParser, a PDBParser object is created, and the structure is retrieved using the get_structure() method.
# A for loop is constructed to iterate over and get the residues in the structure.
# Two varibales are set, residue_1 retieves the x,y,z coordinates of the C-alpha atom of an amino acid residue, using the get_coord() method.
# The variable residue_2 retrieves the x,y,z coordinates of the C-alpha atom of another amino acid residue, this is also done using the get_coord() method.
# The C-alpha distance is calculated using numpy.linalg.norm from the numpy module.
# The numpy.linalg.norm() function is used for calculating different matrix or vector norms.
# The C-alpha distance is calculated and the result is printed out for the user to see.
# The break breaks out of the while loop and another while loop with a try except block is constructed.
# This while loop with the try except block gives the user options of what they would like to do next.
# The user is given the option to calculate C-alpha distance for another PDB file, to return to the Main Menu or exit the program.
# As the input from the user should be an integer associated with the options provided, the except ValueError will detect input that is not an integer from the options provided.
# If user input is not valid then the user will be prompted to input again until there input is valid, this is for input validation.


def function_2():

    input_limit_counter_function_2 = 0

    user_entry = input('Enter the filename of the PDB file you want to work with: ')
    
    while True:

        if input_limit_counter_function_2 == 5:
            print("\nYou have not entered a PDB file within your working directory\n")
            print(
            "You will be redirected to the main page\n"
            + "This will allow you to change your directory as you probably have no PDB files in this directory\n\n")
            time.sleep(3)
            Redirection_to_change_directory()
            break

        if user_entry [-3:] != "pdb":
            input_limit_counter_function_2 +=1
            user_entry = input("This is not a PDB file, please enter a PDB file: ")
            continue

        if not os.path.exists(user_entry):
            input_limit_counter_function_2 += 1
            print("The PDB file %s does not exist in your working directory"%user_entry)
            user_entry = input('Please try again: ')
            continue
        
        
        if user_entry[-3:] == "pdb" and os.path.exists:

            parser = PDBParser()
            structure = parser.get_structure(user_entry, user_entry)

            retrieving_residues = [residues for residues in structure.get_residues()]       

            residue_1 = retrieving_residues[0]["CA"].get_coord()
            residue_2 = retrieving_residues[1]["CA"].get_coord()

            print("\nC-alpha distance = ",numpy.linalg.norm(residue_1-residue_2))

            # Link to where code was adapted from: https://www.biostars.org/p/374180/
            
            break
        
        else:
            user_entry = input('\nYour input is not a PDB file, please try again: ')

    while True:
        try:
            run_again_choice = int(input(
            "\n Enter 1 to calculate C-alpha distance for another PDB file"
            "\n Enter 2 to return to the Main Menu"
            "\n Enter 3 to EXIT the program\n\n"
            ))
            while (run_again_choice < 1) or (run_again_choice > 3):
                print("\nThis option does not exist !\n\nPlease try again and select from the list of options\n")
                run_again_choice = int(input(
                "\n Enter 1 to calculate C-alpha distance for another PDB file"
                "\n Enter 2 to return to the Main Menu"
                "\n Enter 3 to EXIT the program\n\n"
            ))
            break
        except ValueError:
            print("\nYou have not entered a number from the list, please try again\n")
            continue
        else:
            break
    
    
    if run_again_choice == 1:
        function_2()
    
    elif run_again_choice == 2:
        Main_Menu()
    
    elif run_again_choice == 3:
        exit()

# The function function_3() generates plots of B-factor distribution and Predominant Amino Acid Residues.
# A while loop is constructed, inside the while loop there is a try except block.
# Inside the try block the variable plot_selection provides the user with what plots are available for the user to generate.
# The int(input()) means that the user can only enter an integer input, in this case only the integers 1 or 2 for generating plots.
# The inner while loop checks whether the user has inputted an integer less than 1 or greater than 2, if so options are provided again and the user can retry.
# The except block has an except ValueError.
# Therefore if an integer is not entered by the user,the print statement within the except ValueError will be printed and the user will be notified that they need to enter again.
# Outside the while loop there are two if statements, if plot_selection == 1 and if plot_selection == 2.
# Based on whether the user enters 1 or 2 different plots will be generated.
# In both if statements the user is prompted to enter a PDB structure from the Protein Data Bank without the extension .pdb.
# Using the BioPandas package the users entry is retrived.
# Using the matplotlib libary, if the plot_selection == 1 then a Predominant Amino Acid Residues plot is generated.
# If the plot_selection == 2 then a B-factor Distribution plot is generated.
# After the plot has been generated and the user exits the plot and does not want to generate any more plots, the while loop is constructed to give the user options on what to do next.
# The try except block in the while loop operates in a similar way to the previous try except block in the other while loop in function_3().
# Based on what the user enters, the Main_Menu() function, function_3() or exit() to exit the program are called. 

def function_3():
    
    while True:
        try:
            plot_selection = int(input(
            "\nThese are the selection of plots to choose from:\n\n" +
            "Enter 1 for a Predominant Amino Acid Residues Plot\n" +
            "Enter 2 for a B-factor Distribution Plot\n\n"+
            "Enter the number associated with the plot you want to generate\n"
             ))
            while (plot_selection < 1) or (plot_selection > 2):
                invalid_message = int(input("\nThis option does not exist !\n\nPlease try again and select from the list of options\n"))
                plot_selection
        except ValueError:
            print("\nYou have not entered a number from the list, please try again and select from the list of options\n")
            continue
        else:
            break

   
   
    if plot_selection == 1:

        user_entry = input("\nEnter any PDB structure from the Protein Data Bank, don't put the extension .pdb\nFor e.g. for 3eiy.pdb enter 3eiy\n")

        retrieve_pdb_file = PandasPdb().fetch_pdb(user_entry.strip())

        style.use('ggplot')

        retrieve_pdb_file.df['ATOM']['residue_name'].value_counts().plot(kind='bar')
        plt.title('Predominant Amino Acid Residues')
        plt.xlabel('Amino Acid')
        plt.ylabel('Number of ATOM entries')
        plt.show()
        
    if plot_selection == 2:
        
        user_entry = input("\nEnter any PDB structure from the Protein Data Bank, don't put the extension .pdb\nFor e.g. for 3eiy.pdb enter 3eiy\n")
    
        retrieve_pdb_file = PandasPdb().fetch_pdb(user_entry.strip())

        style.use('seaborn-pastel')

        retrieve_pdb_file.df['ATOM']['b_factor'].plot(kind='line')
        plt.title('B-factor Distribution')
        plt.xlabel('Atom Number')
        plt.ylabel('B-factor')
        plt.show()

    while True:
        try:
            run_again_choice = int(input(
            "\n Enter 1 to generate another plot for the same or another PDB file"
            "\n Enter 2 to return to the Main Menu"
            "\n Enter 3 to EXIT the program\n\n"
            ))
            while (run_again_choice < 1) or (run_again_choice > 3):
                print("\nThis option does not exist !\n\nPlease try again and select from the list of options\n")
                run_again_choice = int(input(
                "\n Enter 1 to generate another plot for the same or another PDB file"
                "\n Enter 2 to return to the Main Menu"
                "\n Enter 3 to EXIT the program\n\n"
                ))
            break
        except ValueError:
            print("\nYou have not entered a number from the list, please try again\n")
            continue
        else:
            break

    if run_again_choice == 1:
        function_3()
    elif run_again_choice == 2:
        Main_Menu()
    elif run_again_choice == 3:
        exit()



# The code inside function_4() retireves atomic coordinates x,y,z for an amino acid that is selected by the user.
# Then an output file is generated which can be visualised in PyMOL as PyMOL is opened from Python to allow the user access to the output file/files that the program generates.
# It is called function_4() as it is listed as choice number 4 when running the program.
# The user is asekd to enter a PDB file they want to work on for this function using the input() prompt. 
# Storing all of the three letter codes for the 20 amino acids in a list as later the if conditional will check if the user has entered a three letter code that is in this list.
# The variable input_limit_counter_function_4 is initilised, this variable limits user input.   
# A while loop is used for all the conditional statements relating to how the function operates in relation to user input. 
# If the user enters a three letter code for an amino acid that is present in the amino_acid_codes list then all the atom entries for that amino acid will be shown.
# An output file will be saved with a specific file name relating to the amino acid the user selected.
# A print statement is used to tell the user that an output file has been saved and the name of the output file is shown to the user.
# This allows the user to find the output file and then rename it if they wish to.
# The user is prompted whether they would like to carry the function out again for another amino acid, this is done using the input() prompt.
# The user should enter Y to carry out the function again or N for no.
# If the user does not enter Y or N they will be prompted with the same message that they should enter Y for yes and N for no, this is to validate input.
# If the user enters N, the print statement will show a message telling them that they are being directed to PyMOL.
# The time.sleep(3) from the time module in Python is used to create a three second delay so that the user is aware of what's going to happen and dosen't become overwhelmed with the program.
# The os.system('pymol') from the os module is used as a command line argument which opens PyMOL.
# After the user has finished visualising the output file/files in PyMOL and they exit PyMOL, the input prompt is used to give them a choice of what to do next.
# By using if conditional statements inside a while loop, if the user enters E, sys.exit(0) from the sys module is used to exit the program.
# If the user enters M then the program will return to the Main Menu.
# If the user enters something or than E or M they will be prompted that there input was not an option provided and will be prompted to input until they select an option given by the program.

def function_4():
    
    user_entry = input('Enter the filename of the PDB file you want to work with: ')

    amino_acid_codes  = [
        "ALA","ARG","ASN","ASP",
        "CYS","GLU","GLN","GLY",
        "HIS","ILE","LEU","LYS",
        "MET","PHE","PRO","SER",
        "THR","TRP","TYR","VAL"
        ]
    
    input_limit_counter_function_4 = 0

    while True:

        if input_limit_counter_function_4 == 5:
            print("\nYou have not entered a PDB file within your working directory\n")
            print(
                "You will be redirected to the main page\n"
            + "This will allow you to change your directory as you probably have no PDB files in this directory\n\n"
            )
            time.sleep(3)
            Redirection_to_change_directory()
            Main_Menu()
            break

        elif user_entry [-3:] != "pdb":
            input_limit_counter_function_4 +=1
            user_entry = input("This is not a PDB file, please enter a PDB file: ")
            continue
        
        elif not os.path.exists(user_entry):
            input_limit_counter_function_4 +=1
            print("The PDB file %s does not exist in your directory"%user_entry)
            user_entry = input('Please try again: ')
            continue
        
        elif user_entry[-3:] == "pdb" and os.path.exists:
            PDB_File = open(user_entry)

            user_input = input(
            "\nPlease enter the three letter code of the amino acid (in Capitals) you want to view atomic coordinates for:\n" +
            "for e.g. for GLYCINE enter GLY\n")
            while True:
                if user_input not in amino_acid_codes:
                    user_input = input("That is not a three letter amino acid code, please try again ")
                else:
                    break

            for line in PDB_File:
                if line.startswith("ATOM") and line.__contains__(user_input):
                    print(line.strip())
            
            
            
            print(line.strip(),file=open(user_entry + "_" + user_input + "_outfile.pdb","a"))
                    
            print(
                "\nYour output has been successfully saved as a PDB file named\n"
                + str(user_entry) + "_" + user_input + "_outfile.pdb")
            
            PDB_File.close()
        
        break
        
    run_again = input("\nWould you like to do this for another amino acid?\n" + 
    "Enter Y for yes or N to open PyMOL and visualise your output file/files\n") 

    while True:

        if run_again == "Y":
            function_2()
            break
        if run_again == "N":
            break

        if run_again != "Y" or run_again != "N":
            run_again = input("\nPlease enter Y or N: ")
    
    print("\n\nYou are now being directed to PyMOL where you can view the output file/files you have generated in this function")

    time.sleep(3)
    
    os.system('pymol')

    return_or_exit = input("\nWould you like to return to the main menu or exit the program\n"
        + "Enter E to exit the program or M to return to the main menu\n")

    while True:

        if return_or_exit == "E":
            print("\nThank you for using the program")
            sys.exit(0)
        if return_or_exit == "M":
            print("\nTaking you back to the main menu...\n\n")
            time.sleep(3)
            Main_Menu()
        
        if return_or_exit != "E" or return_or_exit != "M":
            print("\nYou have not selected a option from the options provided! Please try again\n")
            return_or_exit = input("\nWould you like to return to the main menu or exit the program\n"
            + "Enter E to exit the program or M to return to the main menu\n")
    


# The function Main_Menu() contains the code for the Main Menu which gives the user the option to select what functionality they want to make use of.
# Inside the while loop is the try except block which is used for user input and user input validation.
# Inside the try, the segment of code prompts the user for input on which function they want to select based on the number assigned to that function and also the number to select if they wish to exit the program.
# If any number from 1 to 4 is selected the break will exit the while loop and move to the next section of code outside the while loop.
# The expect ValueError is used to only allow user input to be an integer.
# Anything other than a integer entered by the user will be detected and the user will be prompted by the program to enter again.
# The segment of code conssisting of if statements calls functions based on what the user enters.
# So for e.g. if the user has entered 1, function_1() will be called and the code for that will be executed accrodingly.
# In the Main Menu there is also an option where the user enters 5 and the program exits which is one of the if conditionals below.

def Main_Menu():
 
    while True:
        try:
            user_menu = int(input(
                "\n\n Main Menu\n\n"
                " Enter 1 for Structural information of the protein\n" +
                " Enter 2 for C-alpha distance between two residues\n" +
                " Enter 3 for Plotting B-factor Distribution and Predominant Amino Acid Residues\n"
                " Enter 4 for retrieving atom x,y,z coordinates for a selected amino acid and visualisation on PyMOL\n" 
                " Enter 5 to quit and exit the program\n\n"))
            while (user_menu < 1) or (user_menu > 5):
                print("\nThis option does not exist !\n\nPlease try again and select from the list of options\n")
                user_menu = int(input(" Enter 1 for Structural information of the protein\n" +
                " Enter 2 for C-alpha distance between two residues\n" +
                " Enter 3 for Plotting B-factor Distribution and Predominant Amino Acid Residues\n"
                " Enter 4 for retrieving x,y,z coordinates of atoms for a selected amino acid and visualisation on PyMOL\n" 
                " Enter 5 to quit and exit the program\n\n"))
            break
        except ValueError:
            print("You have not entered a number from the list, please try again ")
            continue
        else:
            break


    if user_menu == 5:
        quit()

    if user_menu == 1:
        function_1()
    

    if user_menu == 2:
        function_2()
    
    
    if user_menu == 3:
        function_3()        

    if user_menu == 4:
        function_4()



# The directory_path() function is the first user interactive function of the program.
# After the current working directory is displayed to the user, the while loop handles user input to stay or change the current working directory.
# The if, elif and else conditionals handle input based on what the user inputs.
# If the user enters STAY, the break breaks out of the while loop and the Main Menu is displayed.
# If the user enters CHANGE, an input prompt will allow the user to enter the path of the directory they want to become their new working directory.
# A try except block is used in the elif statement to so that the path the user enters is a path that exists. 
# If the path is not a path that exists, the OSError in the except block will tell the user that the path does not exist.
# The user will then be prompted to enter the path again.
# In the try block, if the path that the user enters does exist, the current working directory will be changed.
# The print statement in the try block tells the user that the working directory has changed and prints the path of the new working directory.
  
def directory_path():

    while(True):
        change_wd = input("\n\nEnter STAY to remain in the current working directory or CHANGE to change the current working directory\n")
        if change_wd == "STAY":
            break
        elif change_wd == "CHANGE":
            dir_change = input("\nEnter the path of the directory you would like to become your new working directory\n") 
            try:
                os.chdir(dir_change)
                print("Working directory has changed to:\n\n" + dir_change)
                break
            except OSError:
                print("\nThat is not a path that exists\n")
                dir_change
        else:
            print("\nNot a valid option, please try again\n")
            change_wd



# Two print statements, one to welcome the user to the program and the other above that to format the print statement moving it two lines below. 
# This is for better visualisation of the welcome print statement in the command line.
# The ANSI escape sequecnes are used on the welcome print statement to make it bold so the print statement acts as the header of the program.
# The ANSI escape sequences are the segments \033[1m on either end of the second print statement below.

print("\n\n")
print("\033[1m                              Welcome to the PDB File Parser Program                         \033[1m")



# The current working directory is obtained by using the os module, the variable working_directory_path is set to os.getcwd() which gets the current working directory.
# The print statement prints the string shown below and the variable working_directory_path. This outputs the current working directory path.

working_directory_path = os.getcwd()

print("\n\nYour current working directory is: " + working_directory_path)


# Both directory_path and Main_Menu() functions are called.
# This moves into the main frame of the program, allowing the user to either stay or change the current working directory and then select what functionality they want to make use of.
 
directory_path()
Main_Menu()
