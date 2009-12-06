#!/usr/bin/env python
#
#  ARFF.py
#  
#
#  Created by Chris Johnson on 11/10/09.
#  Copyright (c) 2009 The University of Tulsa. All rights reserved.
#

import csv
import re
from itertools import chain,izip as zip

def transpose(matrix):
	"""Transpose a matrix. Assumes an m x n matrix."""
	return [list(m) for m in zip(*matrix)]

def genAttr(attribute):
	if attribute[0] in "'\"":
		raise Exception("Complex attribute names not implemented")

	name, attr = attribute.split(None,1)
	if re.match("numeric",attr,re.I):
		return (name, Numeric())
	if re.match("string",attr,re.I):
		return (name, String())
	if attr[0]+attr[-1] == "{}":
		return (name, Nominal(*(attr[1:-1].split(','))))
	raise Exception("Data attributes not implemented")

# ARFF Attribute types
class Attribute(object):
	"""Base attribute class. Will only match "?" without throwing an error."""
	def match(self, value):
		return value=="?" or re.match(self.pattern, value) is not None
	def __str__(self):
		raise NotImplemented("abstract")

class Numeric(Attribute):
	"""Numeric attribute class. Matches integer or floating representations."""
	pattern = re.compile("(\d+\.?|(\d+)?\.\d+)((e|E)(\-|\+)?\d+)?$")
	
	def __str__(self):
		return "numeric"
	
class Nominal(Attribute):
	"""Nominal attribute class. Matches among a list of strings."""
	def __init__(self, *args):
		self.valid = args
		
	def __str__(self):return "{" + ','.join(self.valid) + "}"
	def match(self, value):return value=="?" or value in self.valid

class String(Attribute):
	def __init__(self, q="'"):
		self.pattern = re.compile(r"%s(([^\\%s]|\\[^%s])*\\%s)*([^\\%s]|\\[^%s])*%s$|[^%s\s]\S*$" % (q,q,q,q,q,q,q,q))
		
	def __str__(self):
		return "string"
		
class Date(Attribute):
	def __str__(self):
		return "Unimplemented"

class ARFF:
	p_rel = re.compile("@relation\s+",re.I)
	p_attr = re.compile("@attribute\s+", re.I)
	p_data = re.compile("@data\s*$", re.I)

	def __init__(self, name = "", attributes = [], data = [[]], comments = []):
		self.name = name
		self.attributes = attributes
		self.data = data
		self.comment = comments

	def __str__(self):
		output = []
		output+=("%% %s" % comment for comment in self.comment)
		output.append("@relation %s" % self.name)
		output+=("@attribute %s %s" % attr for attr in self.attributes)
		output.append("@data")
		output+=(','.join(row) for row in self.data)
		return '\n'.join(output)
		
	def validate(self):
		if not String().match(self.name):
			raise ValueError("ARFF relation name must be defined.")
		if not all([isinstance(attr, tuple) for attr in self.attributes]):
			raise TypeError("Attributes must be tuples.")
		elif not all([isinstance(attr[0], str) for attr in self.attributes]):
			raise TypeError("First element of each attribute must be a string.")
		elif not all([isinstance(attr[1], Attribute) for attr in self.attributes]):
			raise TypeError("Second element of each attribute must be a set.")

		dataT = transpose(self.data)

		if len(dataT) != len(self.attributes):
			raise ValueError("Data rows must have same length as attributes.")
		elif not all([attr.match(val) for (_,attr), values in zip(self.attributes, dataT) for val in values]):
			raise ValueError("Data rows may not have unexpected values.")
		if not all(isinstance(comment,str) for comment in self.comment):
			raise ValueError("Comments must be strings.")

	def writeARFF(self,file):file.write(self)

	def writePARF(self,file):
		output = []
		output+=("%% %s" % comment for comment in self.comment)
		output.append("@relation %s" % self.name)
		output+=("@attribute %s %s" % attr for attr in self.attributes)
		output.append("@data")
		for row in self.data:
			output+=(','.join(partial_row) for partial_row in self.splitrow("?" if x=="NA" else x for x in row))
		output.append("")
		file.write('\n'.join(output))

	def splitrow(self, row):
		if len(row) < 500:return [row]
		return [row[:500] + ["&"]] + self.splitrow(row[500:])
	def splitrow(self,row):#No longer recurses and now capable of using a generator
		row=iter(row).next
		r=[]
		try:
			while 1:r+=chain((row() for a in xrange(500)),"&")
		except StopIteration:return r

	def parseARFF(self,file):
		for line in file:
			line = line.strip()
			if not line:continue

			if line[0]=='%':
				self.comment.append(line[1:])
			elif ARFF.p_rel.match(line):
				self.name = ARFF.p_rel.sub("",line)
			elif ARFF.p_attr.match(line):
				self.attributes.append(genAttr(ARFF.p_attr.sub("",line)))
			elif ARFF.p_data.match(line):
				break

		reader = csv.reader(file)
		self.data = [row for row in reader]
		self.validate()

class BioData:
	def __init__(self, subjects = [], snps = [], genotypes = [], phenotypes = []):
		self.subjects = subjects
		self.snps = snps
		self.genotypes = genotypes
		self.phenotypes = phenotypes

if __name__ == '__main__':
	a = ARFF()
	a.parseARFF(open("/Users/administrator/Downloads/dat2-final.arff",'r'))
