#
#  ARFF.py
#  
#
#  Created by Chris Johnson on 11/10/09.
#  Copyright (c) 2009 The University of Tulsa. All rights reserved.
#

import re

def transpose(matrix):
	"""Transpose a matrix. Assumes an m x n matrix."""
	return [list(m) for m in zip(*matrix)]

# ARFF Attribute types
class Attribute(object):
	"""Base attribute class. Will only match "?" without throwing an error."""
	def match(self, value):
		if value == "?":
			return True
		return re.match(self.pattern, value) is not None

	def __str__(self):
		abstract

class Numeric(Attribute):
	"""Numeric attribute class. Matches integer or floating representations."""
	pattern = re.compile("(\d+\.?|(\d+)?\.\d+)((e|E)(\-|\+)?\d+)?$")
	
	def __str__(self):
		return "numeric"
	
class Nominal(Attribute):
	"""Nominal attribute class. Matches among a list of strings."""
	def __init__(self, *args):
		self.valid = args
		
	def __str__(self):
		return "{" + ','.join(self.valid) + "}"
		
	def match(self, value):
		if value == "?":
			return True
		return value in self.valid
		
class String(Attribute):
	def __init__(self, q="'"):
		self.pattern = re.compile(r"%s(([^\\%s]|\\[^%s])*\\%s)*([^\\%s]|\\[^%s])*%s$|[^%s\s]\S*$" % (q,q,q,q,q,q,q,q))
		
	def __str__(self):
		return "string"
		
class Date(Attribute):
	def __str__(self):
		return "Unimplemented"

class ARFF:
	def __init__(self, name = "", attributes = [], data = [[]], comments = []):
		self.name = name
		self.attributes = attributes
		self.data = data
		self.comment = comments

	def __str__(self):
		output = []
		output.extend(["%% %s" % comment for comment in self.comments])
		output.append("@relation %s" % self.name)
		output.extend(["@attribute %s %s" % attr for attr in self.attributes])
		output.append("@data")
		output.extend([','.join(row) for row in data])

		return '\n'.join(output)
		
	def validate(self):
		if self.name == "":
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

		if not all([isinstance(comment,str) for comment in comments]):
			raise ValueError("Comments must be strings.")

	def writeARFF(self,file):
		file.write(self)

def BioData:
	def __init__(self, subjects = [], snps = [], genotypes = [], phenotypes = []):
		self.subjects = subjects
		self.snps = snps
		self.genotypes = genotypes
		self.phenotypes = phenotypes

	
