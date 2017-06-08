from Bio import SwissProt, ExPASy, SeqIO
import urllib
import urllib2
import argparse
import re


class COGs:
	"""COGs class find the orthologous clusters connected to a protein."""

	def __init__(self, uprotID):
		self.uprotID = uprotID

	def prot_kegg(self):
		"""Transfer uniprot ID to KEGG."""
		uprot = self.uprotID
		kegg = urllib.urlopen("http://rest.kegg.jp/conv/genes/uniprot:"+uprot)
		string = kegg.read()
		keggID = string[string.find("up:")+10:]
		return keggID

	def find_COG(self):
		"""Find COG output for KEGG ID."""
		query = self.prot_kegg()
		url_open = urllib.urlopen("http://rest.genome.jp/oc/?"+query)
		records = url_open.read()
		return records

	def __str__(self):
		"""Prints out keggID and uniprotID."""
		print("UniprotID: %s KeggID: %s" % (self.uprotID, self.prot_kegg()))

	def find_COG2(self):
		"""Find records from uniprotIDs without use of keggIDs."""
		handle = ExPASy.get_sprot_raw(self.uprotID)
		record = SwissProt.read(handle)
		query = record.gene_name.strip("Name""="";")
		url_open = urllib.urlopen("http://rest.genome.jp/oc/?"+query)
		return url_open.read()

	def make_COGlib(self):
		"""Creates a nested list for COG record of a gene."""
		record = self.find_COG()
		match = record.split("\n")		# Make a list with lines in record as elements.
		list_ = []
		for line in match:		# Separate line in categories and make array of record.
			match2 = line.split("\t")
			match3 = re.findall(r'^.(...:.*)\.*$', line)
			list_.append(match2[0:7]+match3)
		return list_

	def cog_ids(self):
		"""Return list with cog ids KEGGid of search-gene and
		length for a gene cluster."""
		dic = self.make_COGlib()
		list_2 = []
		try:
			j = dic[2][0]		# Find OCid.
			list_2.append(j)
			for i in dic:
				if i[0] != j:
					list_2.append(i[0].replace('# ', ''))
			list_2.append("nr genes in cluster: %i" % (len(dic)-2))
			del list_2[2:4]		# Remove unnecessary elements.
		except:
			j = 'no record found for %s' % (self.uprotID)
			list_2.append(j)
		return list_2


class Disease:
	"""Disease class takes a disease as string and do a uniprot query for
	proteins connected to the disease then create a dictionary of gene ID
	keys and swissprot record values."""

	def __init__(self, dis):
		self.dis = dis

	def IDs(self):
		"""Return a list of IDs from query."""
		query_dict = {'query': 'disease'+" "+self.dis, 'sort': 'score',
												'columns': 'id', 'format': 'list'}
		encode_str = urllib.urlencode(query_dict)
		url_open = urllib.urlopen("http://www.uniprot.org/uniprot/?" + encode_str)
		IDs_object = url_open
		ID_list = []
		for i in IDs_object:
			i = i.strip('\n')
			ID_list.append(i)
		return ID_list

	def records(self):
		"""Return a dictionary of ID and swissprot records from query."""
		record_dict = {}
		except_ids = []
		for i in self.IDs():
			try:
				handle = ExPASy.get_sprot_raw(i)
				record_dict[i] = SwissProt.read(handle)
			except HTTPError, AssertionError:
				print("there was a problem finding uniprotID {} \n\
					try Records_fromfile-method".format(i))
		return record_dict

	def records_fromfile(self):
		"""For use with poor network abilities, takes record from Uniprot file."""
		records_ofile = {}
		dictionary = SeqIO.index("uniprot_sprot.dat", "swiss")
		for i in self.IDs():
			records_ofile[i] = dictionary[i]
		return records_ofile

	# This is a bit slow.
	def find_cogs(self):
		"""Find cogs for all proteins associated to a disease."""
		ids = self.IDs()
		length = len(ids)
		cogs = []
		for i in range(length):
			cog = COGs(ids[i])
			cogs.append(cog.cog_ids())
		return cogs

	# Unfinished method, no output.
	def common_cogs(self):
		"""Return uniprot IDs that are part of same COG."""

		ids = self.IDs()
		records = []
		for i in ids:
			cog = COG(ids[i])
			records.append(cog.make_COGlib())
		for i[0] in records:
			if i[0] in records:		# Find duplicate OC.
				KeggID.append(i[1])		# Append keggid name in list.
		return KeggID


class Comparison:
	"""Comparison class compares protein associated uniprotIDs of two diseases."""

	def __init__(self, Disease_1, Disease_2):
		self.Disease_1 = Disease_1
		self.Disease_2 = Disease_2

	def same_protIDs(self):
		"""Return a list of similar proteinIDs if they exist."""
		compared = []
		list1 = Disease(self.Disease_1).IDs()
		list2 = Disease(self.Disease_2).IDs()
		if len(list1) > len(list2):
			longest = list1
			shorterst = list2
		else:
			longest = list2
			shortest = list1
			for i in longest:
				if i in shortest:
					compared.append(i)
		if len(compared) < 1:
			return "no similar proteinIDs"
		else:
			return compared

	# Unfinished method, no output.
	def same_cogsID(self):
		"""Compare cogs for associated proteins of two diseases."""
		d1 = Disease(self.Disease_1)
		d2 = Disease(self.Disease_2)
		d1_cogs = d1.find_cogs()
		d2_cogs = d2.find_cogs()
		compared = []
		for i in d1_cogs[0]:
			if i in d2_cogs[0]:
				compared.append(i)
		if len(compared[0]) < 1:
			return "no similar cogs"
		else:
			return compared[0]


def main():
	"""Function for command line interface."""
	parser = argparse.ArgumentParser()
	parser.add_argument("string", nargs='+', help='give disease to retrieve \n\
		uniprot IDs of associated proteins')
	parser.add_argument("-c", "--cogfind", action='store_true', help='display a \n\
		the cog ids for the associated proteins of a disease')
	parser.add_argument("-s", "--simprot", action='store_true', help='display \n\
		uniprotIDs that are common in two diseases')
	args = parser.parse_args()
	input_disease = args.string

	if args.cogfind:
		if len(input_disease) == 1 :
			answer = Disease(input_disease[0])
			print(answer.find_cogs())
		else:
			print("One ID at a time please, e.g: python Record_dict.py disease ")

	elif args.simprot:
		if len(input_disease) < 2:
			print("need two diseases e.g: python Record_dict.py -s disease1 disease2")
		else:
			in1 = input_disease[0]
			in2 = input_disease[1]
			print(Comparison(in1, in2).same_protIDs())

	else:
		if len(input_disease) == 1 :
			answer = Disease(input_disease[0])
			print(answer.IDs())
		else:
			print("One ID at a time please, e.g: python Record_dict.py disease ")


if __name__ == '__main__':
	main()
