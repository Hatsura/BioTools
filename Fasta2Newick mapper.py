'''
	'FASTA header into Newick' Mapper
	release date: 17-jan-2018
	by Paulo 
	GitHub: @Hatsura
	
	version: 1.1
	last updated: 23-01-2018
'''
# -- CONFIG--
# This separator character will be used to separate description from specific name.
# ie.: homeobox protein - Apis melifera
# IMPORTANT NOTE: you can't use any of these following characters: , ( ) ;
sep = '@' 

# -- END of CONFIG --

# Aesthetics: adding spaces around separator
sep = ' '+sep+' '

#Foolproofing
if sep in [',',':',';',')','(']:
	sep = ' '

from tkinter import Tk,filedialog #Stuff for file open/saving. -- https://pythonspot.com/en/tk-file-dialogs/

def parse_header(header):
	'''
	Returns the Acession ID, protein description and organism name
	'''
	# Removing newlines
	header = header.strip('\n').strip()
	
	# First, get the organism name
	temp_header = header.split('[')				# ID desc + Organism name
	organism_name = temp_header[1].replace(']','').strip()	# Organism name, removing '\n'
	# Removing organism name from splitted header's list.
	# Acession and Description remains.
	temp_header.remove(temp_header[1]) 

	# Second, get the ID.
	# Separating by spaces, acession will always be the first...in current NCBI format.
	acession = temp_header[0].split(' ')[0].replace('>','').strip()

	# The remaining text is the description
	# Cleaning process: removal of acession, removal of surrounding spaces and removing
	# the NEWICK-illegal chars (':', '(', ')'' )for '--'.
	description = temp_header[0].replace(acession,'').strip(' ').replace(':','--').replace('(','--').replace(')','--')

	return (acession, description, organism_name)

if __name__ == '__main__':
	root = Tk() 	# Main window (root)
	root.withdraw() # Hiding 'root' window
	My_file = filedialog.askopenfilename(title='Select the FASTA file containing the sequence headers.',
	                                      filetypes = ( ("FASTA files \(*.fasta\)","*.fasta"),
	                                                    ("All files","*.*")
	                                                  )
	                                      )
	# Parsing file
	with open(My_file,'r') as f:
		headers = f.readlines()

	# Filtering and leaving only the headers
	# Also, removing the '>'
	headers = [x.replace('>','') for x in headers if '>' in x]

	parsed_headers = []
	for h in headers:
		parsed_headers.append(parse_header(h))

	My_file = filedialog.askopenfilename(title='Select your NEWICK format file.',
	                                      filetypes = ( ("NEWICK files \(*.nwk\)","*.nwk"),
	                                                    ("All files","*.*")
	                                                  )
	                                      )

	with open(My_file,'r') as f2:
		phylogram = f2.readlines()
	
	#Converting into a string
	phylogram = ''.join(phylogram)

	# Remapping IDs -> description [organism name]
	for k in parsed_headers:
		name = k[1] + sep + k[2]

		# In case other protein has same 
		if name in phylogram:
			# Adding a cumulative '*' for same node proteins.
			name += '*'	
		phylogram = phylogram.replace(k[0], name)

	result_file = filedialog.asksaveasfilename( title='Saving converted NEWICK file into...',
	                                                defaultextension='.nwk',
	                                                filetypes= [('*.nwk','.nwk')]
												)

	with open(result_file, 'w') as f3:
		f3.write(phylogram)