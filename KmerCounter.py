#!/usr/bin python3
# *-* coding: utf-8 *-*
''' versão: 2.0.0
 por: Paulo
 GitHub: @Hatsura
 
 - Script feito quando ainda estava na graduação e não manjava muito de bioinfo ou de programas melhores para fazer isso (Jellyfish).
 - Conta os kmers numa sequência, retorna suas frequências relativas achadas na sequeência. 
 - A configuração do tamanho dos k-mers a ser analisado tem de ser fornecido manualmente no script (buscar por !CONF).
 - Por padrão, analisa os 3 frames de leitura (janela deslizante com step 1). Dá para modificar fornecend o parâmetro 'Step' (buscar !CONF2)
 - Aceita single e multi fasta. Caso detecte várias sequências, fornece opção de salvar em arquivo único ou múltiplos arquivos.
 - A única parte com interface gráfica são as janelas de carregar e salvar arquivo. :)
 - Não tem linha de comando. Basta configurar e depois executar o script.
 - Boa parte dos comentários do código estão em inglês. Por motivos de porque sim.
'''

from tkinter import filedialog,Tk #Stuff for file open/saving. -- SRC: https://pythonspot.com/en/tk-file-dialogs/
import errno, random, itertools, sys
from os.path import exists, dirname
from os import makedirs
import pandas as pd

'''
    IUPAC nucleotide Ambiguity Codes
------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
A            A           Adenine
C            C           Cytosine
G            G           Guanine
T            T           Thymine
U            U           Uracil
M          A or C
R          A or G
W          A or T
S          C or G
Y          C or T
K          G or T
V        A or C or G
H        A or C or T
D        A or G or T
B        C or G or T
X      G or A or T or C
N      G or A or T or C
'''
def Read_Seq( Seq,
              Seq_Name = 'default_sequence_name',
              #Reversion_Check = False,
              Oligo_sizes = [1],
              Only_ATCG = True,
              Step = 1,
              RSCU = False):
    '''
    Seq = Input sequence
    Oligo_sizes = list containing desired sizes of oligonucleotides \
                    to look for in the following. \
                    By default, checks for oligos ranging from size 1 to 4.
    Only_ATCG = Boolean flag on wether the user wants to also analyze non-ATCG\
                nucleotides (good for poorly-made sequences, which shouldn't be
                usefull anyway)
    Reversion_check = ?
    Step = Step of the sliding window. Default = 1.
    '''

    # Exporting variables to log file.
    global cfg_oligos, cfg_step
    cfg_oligos = Oligo_sizes
    cfg_step = Step

    if Step not in [1,2,3]:
        # Not between 1~3 but is positive
        if Step >= 0:
            print("Error: Step size must be within 1~3 range.")
            print("Use negative numbers to override it.")
        elif Step < 0:
            Step = -Step # Turning Step size positive. FIXME: Should raise an error here instead?
        
    # Checking the k-mer sizes to be checked. Also, removing duplicates.
    if 0 in Oligo_sizes:
        Oligo_sizes = list(set( [x for x in Oligo_sizes if x != 0] ))
        print('Warning: Ignoring \'0\' (zero) as oligonucleotide size to be analized.')

    # -- Declaring Variables --
    DNA_Data = {}   #Dict that will hold whole genome data. Will be converted
                    #to pd.DataFrame
   
    #Sequence Information-related
    Sequence_Size = len(Seq)# Genome Size!
    Valid_Nucleotides = []    # Will store the valid nucleotides 
                            # to populate dictionary and to be look for
    Oligo_list = [] # Holds the oligonucleotides that will be iteratively 
                    # added to the 'DNA_Data' dictionary.
    Kmer_count = 0  #Stores the count of k-sized kmers found for each oligo size defined by user.

    # User can chose between analyzing also non-ATCG nucleotides or -only- ATCG 
    # nucleotides (default)
    if Only_ATCG == True:   
        Valid_Nucleotides = ['A','T','C','G']

    else: #User wants to also analyze non-ATCG nucleotides. Go figure...
        for n in ['A','T','C','G','Y','R','W','S','K','M','D','V','H','B','N','X','U']:
            if Seq.count(n) > 0:
                Valid_Nucleotides.append(n)
            
        # Checks and tell if/how many non-ATCG nucleotides were found more than just ATCG nucleotides 
        # in order to warn the user about it.
        #FIXME: Não tá funcionando!!
        print('\n')
        print("[WARNING] Non-ATCG nucleotides ("+str(len(Valid_Nucleotides))+") were found in this sequence: "  \
             + ' '.join( list( set(Valid_Nucleotides) - set(list(['A','T','C','G']) ) ) ) 
            )

    # Sanity Check. Pretty impossible a DNA with 3 or less nucleotidic bases.
    if len(Valid_Nucleotides) < 4:
        print('[WARNING]: Your sequence has less than 4 types of nucleotides: ' \
               + ' '.join(Valid_Nucleotides) )
    
    # Converting a list into this string so it can be used below 
    # in the 'itertools.product' function.
    Valid_Nucleotides = ''.join(Valid_Nucleotides)

    #Populating 'DNA_Data' dictionary...
    for i in sorted(Oligo_sizes): 
        
        # Populates the dictionary only with the oligonucleotides 
        # provided in the 'Oligo_sizes' list, located right below.
        Oligo_list = list(itertools.product( Valid_Nucleotides,
                                             repeat = i)
                                             )
        
        # Declaring all oligos into a dictionary to prevent errors.
        # Also, it's faster than checking every single iteration if it already exists. (wew, i think so... :D)
        for j in Oligo_list: 
            DNA_Data[ ''.join(j) ] = {'Count': 0}

    # Data Gathering...START!
    # -- Sequence Analysis -- 
    # Time to sweep the given genomic sequence looking for nucleotides!
    # For each nucleotide at index 'i' in the genomic sequence...
    for i in range( 0, Sequence_Size , Step):
        # Then look for 'kmer_size' of 'i+kmer_size' oligo size.
        for kmer_size in Oligo_sizes: #Increment the count of found oligonucleotide of size 'j'
            if i + kmer_size > Sequence_Size:
                continue    #Can't extrapolate the string! Try next oligo size.

            # Incrementing k-mer count
            else:
                # Sequence present in current sliding window
                window_seq = Seq[i: i+kmer_size]
                try:
                    DNA_Data[ window_seq ]['Count'] += 1
                    Kmer_count += 1
                except:
                    continue
    
    # User wants the RSCU of this data.
    # Returning a dataframe of the RSCU
    if RSCU == True:
        return pd.DataFrame.from_dict( GetRSCU( in_dict = DNA_Data,
                                                seqname = Seq_Name),
                                       orient='index')
    
    else:
        # Calculating the relative frequencies and adding it into a separate column!
        # 
        # Formula for obtaining the total amount of k-mers of size 'k' found: 
        #                    Frel = C / (S - (k -1))
        #     Where:
        #        C = count (frequency) of a given nucleotide/aminoacid
        #        S = sequence size
        #        k = k-mer length (size)
        for i in DNA_Data:
            if Step == 1:
                DNA_Data[i][Seq_Name] = round( DNA_Data[i]['Count'] / float( Sequence_Size - (len(i) - 1) ), 6)
            # In case user used a Step greater than 1...
            # FIXME: Wouldnt be easier/universally secure to just use this instead?
            else:
                DNA_Data[i][Seq_Name] = round(( DNA_Data[i]['Count'] / Kmer_count), 6)
        
        # Converting: Dict -> DF

        DNA_Data = pd.DataFrame.from_dict(DNA_Data, orient='index')
        
        # You don't need the counts in the relative frenquencies
        # TODO: add a bool for this? Make it user's-choice?
        DNA_Data = DNA_Data.drop('Count', axis=1)

        # Returning a dataframe with the relative frequencies of k-mer counts.
        return DNA_Data
    
# -----------------------------------------
# -- Removes junk words from the header  --
# -----------------------------------------
# TODO: Make the 'remove' list be loaded from the cfg.ini file as a global var.
# TODO2: Auto-detect if it's a Protein or Genomic sequence.
def clean_header(header):
    # Filtering, removing AID.
    
    # FIXME: What's this?
    out_header = ' '.join(header.split(' ')[1:])
    
    # Don't change anything if there are no spaces. ie.: a simple header (>cas9)
    # Cheapskate: if there are no spaces, there is just 1 text column at header
    if out_header.count(' ') <= 1:
        return header.replace('>','')
    
    # Else, clean it up.
    # Remove ugly terms. ./Absolutely_Barbaric.jpg
    removable = [
            'complete genome',
            ',',
            'genomic sequence',
            '>'
            ]
    
    # Search 'n destroy!
    for i in removable:
        out_header = out_header.replace(i,'')

    return out_header

# ---------------------------------------------------------------------------
# > Helper function:													   --
# 	>> Reads a fasta file, extracts all pairs of 'header + sequence' found --
# ---------------------------------------------------------------------------
def fasta_parse(fasta_path):
    '''
    Returns: a list of lists. 
    Format of each element: [Header, Sequence]
    '''
    # Declaring
    header = ''
    sequence = ''
    Sequences = []
    
    purge = lambda x: x.replace('\n','').replace('\r','')
    
    with open(fasta_path, 'r') as f:
        fastafile = f.readlines()

    for line in fastafile:
        
        if '>' in line:
            # Reached a new header. Save previous, start new parsing/storing.
            if header != '':
                Sequences.append( [
                                purge(header),  # Seq's Header
                                purge(sequence) # Seq's Sequence
                                ]
                            )
                
            # This line is the header        
            header = line
            sequence = '' #Resetting sequence string.
    
        # Tieing up the sequence in a string.
        else:
            sequence += line  
            
    else: #End of loop's else
        Sequences.append( [
                    purge(header),  # Seq's Header
                    purge(sequence) # Seq's Sequence
                    ]
                )
        
    return Sequences

        
#%%==============================================================================
#   -- Performing Genomic Analysis --
#
#   With file & settings decided, time to take action here!
#==============================================================================
def run_analysis(Sequence, Operation, Seq_header):

    data = Read_Seq(Sequence,
    					Oligo_sizes = [1,2,3,4], # !CONF - Tamanhos de kmer a serem analisados. 
    					#Step = 3, # !CONF2  (step=3) == MODO DANIEL -- Optional parameter. Default: 1.
    					Seq_Name = Seq_header,
    					Only_ATCG = True,     #User wants to check only ATCG nucleotides.
					)    
    return data

# =======================
# -- Saves the results --
# =======================
def SaveFile(x, fname = '', silent = False):
    # Prompting the path for the saved data.
    result_file = filedialog.asksaveasfilename(title = 'Saving file into...', 
                                               defaultextension = '.tsv',
                                               initialfile = fname,
                                               filetypes = [('*.tsv','.tsv')]
                                               )
    if result_file:
        # Saving into a '*.tsv' file.
        x.to_csv(path_or_buf = result_file, sep = '\t')
        if not silent:
            print('\t-- Results were saved at:', result_file)

        # Preveting overwriting when saving multiple individual output files.
        # Appending to list
        global cfg_savedfiles
        cfg_savedfiles.append(result_file)
        
    else:
        print('\t-- Results were NOT saved.')    

    return

# ============================================
# -- Prompts a file browser for input files -- 
# ============================================
def GetFiles():
    files = filedialog.askopenfilenames(title='Select file(s) to be analyzed',
                                          filetypes = ( ("FASTA files \(*.fasta\)","*.fasta"),
                                                        ("Text files \(*.txt\)","*.txt"),
                                                        ("All files","*.*")
                                                      )
                                          )

    if files == '':
        print('Error: You must chose a file to work with!')
        sys.exit('')

    return files

# =======================================================================
# -- Should output file be saved separetely or in a single file?       --
# =======================================================================
def GetSaveEach():
    choice = None
    while choice == None:
        choice = input('Multiple input files were detected.\nDo you want to save the results of each input file individually? [y/N] \n> ')
        if choice == '' or choice.lower() == 'n' or choice.lower() == 'no':
            choice = False
        elif choice == choice.lower() == 'y' or choice.lower() == 'yes':
            choice = True

    return choice

#==============================================================================
# -------------------------
#   Main Script Settings
# -------------------------
#==============================================================================

# LOGfile variables
cfg_savedfiles = []

#Initializing GUI
root = Tk()     # Janela 'root', 'root de todos os males'.
root.withdraw() # Esconde a janela 'root'

#  --------------------------
#        Main Script
#  --------------------------
# TODO/FIXME: Send all configurations to a single class called 'configs'? :)
My_data = pd.DataFrame() # This dataframe will hold the resulting data.
My_files = GetFiles()
if len(My_files) > 1:
    save_each = GetSaveEach() #1 output file per input file -or- 1 single output?
else:
    save_each = False # Save everything in a single file.
print('    -- Analyzing',len(My_files),'file(s) --')

#   -- Iterate trough files, extracting k-mer frequencies.
for i in range( len(My_files) ):

    # This will return a list with all stored sequences in the file.
    My_FASTA = fasta_parse(My_files[i])

    # Verbose
    print(str(len(My_FASTA)),'sequences were found in', My_files[i].split('/')[::-1][0])

    # After parsing every sequence in the multifasta, analyze each sequence.
    for header, sequence in My_FASTA:

        # Send each found sequence in this file to analysis.
        temp_data = run_analysis( Sequence = sequence.upper(), #Sequences must be converted to uppercase else it won't be parsed properly.
                                  Operation = 1,
                                  Seq_header = clean_header(header)
                                  )

        # Concatenate all partial data into the main dataframe.
        My_data = pd.concat( [My_data, temp_data],
                             axis = 1)
    
        temp_data = None #Freeing memory

    # Adding values found to the results dataset.
    # Saving now instead of just in the end.
    if save_each == True:
        # Suggested output name. Its based in the input file.
        out_name = My_files[i].split('/')[::-1][0].split('.')[0]

        # Verbose
        #print('Finished analyzing', My_files[i])
		
        # File's Results. Transposing the dataframe so it is in a 'ready2work' format.
        My_data = My_data.transpose()
        
        SaveFile(My_data, fname = out_name, 
                 silent = True)

        #resetting variables/dataframes.
        My_data = pd.DataFrame()
        
# User wants to save all results in a single file.
if save_each == False:
    # Final Results. Transposing the dataframe so it is 'ready2work'.
    My_data = My_data.transpose()

    # Final part: save the data!
    SaveFile(My_data)

logfile = filedialog.asksaveasfilename(title = 'Saving logfile into...', 
                                               defaultextension = '.txt',
                                               initialfile = 'LOG_',
                                               filetypes = [('*.txt','.txt')]
                                               )
if logfile:
    with open(logfile,'w') as f:
        f.write('K-mer sizes analyzed: ' + str(cfg_oligos) + '\n')
        f.write('Sliding window step: ' + str(cfg_step) + '\n' )
        f.write('The following input files were analyzed:' + '\n')
        for i,fin_path in enumerate(My_files,1):
            f.write('\t'+str(i)+'- '+fin_path  + '\n')
        f.write('The files were saved individually? '  + str(save_each) + '\n' )
        f.write('Output file(s): \n')
        for i,fout_path in enumerate(cfg_savedfiles,1):
            f.write('\t'+str(i)+'- '+fout_path + '\n' )