'''
	'FASTA header into Newick' Mapper
	by Paulo 
	GitHub: @Hatsura
	
	version: 1.0
	release date: 17-jan-2018
'''
from tkinter import Tk,filedialog #Stuff for file open/saving. -- https://pythonspot.com/en/tk-file-dialogs/

def parse_header(header):
	'''
	Returns the Acession ID, protein description and organism name
	'''
	# Removing newlines
	header = header.strip('\n').strip()
	
	# First, get the organism name
	temp_header = header.split('[')				# ID desc + Organism name
	organism_name = temp_header[1].replace(']','').strip()	# Organism name
	temp_header.remove(temp_header[1]) #Removing from list

	# Second, get the ID
	acession = temp_header[0].split(' ')[0].replace('>','').strip()

	# The remaining text is the description
	description = temp_header[0].replace(acession,'').replace('>','').strip(' ').replace(':','--')

	return (acession, description, organism_name)

if __name__ == '__main__':
	root = Tk() #Janela 'root', 'root de todos os males'.
	root.withdraw() #Esconde a janela 'root'
	My_file = filedialog.askopenfilename(title='Selecione o seu arquivo FASTA com os cabeçalhos.',
	                                      filetypes = ( ("FASTA files \(*.fasta\)","*.fasta"),
	                                                    ("All files","*.*")
	                                                  )
	                                      )
	# Parsing file
	with open(My_file,'r') as f:
		headers = f.readlines()

	# Filtrando e deixando apenas os cabeçalhos
	headers = [x for x in headers if '>' in x]

	parsed_headers = []
	for h in headers:
		parsed_headers.append(parse_header(h))

	My_file = filedialog.askopenfilename(title='Selecione o seu arquivo NEWICK da árvore.',
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
		name = k[1] + ' ' + k[2]

		# In case other protein has same 
		if name in phylogram:
			name += '*'
		phylogram = phylogram.replace(k[0], name)

	result_file = filedialog.asksaveasfilename( title='Saving converted NEWICK file into...',
	                                                defaultextension='.nwk',
	                                                filetypes= [('*.nwk','.nwk')]
												)

	print(phylogram)
	with open(result_file, 'w') as f3:
		f3.write(phylogram)
