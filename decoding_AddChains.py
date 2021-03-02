#!/usr/bin/env python3
import json
import click
import sys

def substring(s : str, start : str, end : str): return s[s.find(start)+len(start):s.rfind(end)]

def rewrite_line(line : str):

	SPECIAL_SQR = (line[:2] == '2*')

	if '(' in line and ')' in line:

		string = substring(line, '(', ')')
		second = line[len(string) + 2 + 2*SPECIAL_SQR:]
		tmp = string
		#string = ('fp_sqr(' * SPECIAL_SQR) + rewrite_line(string) + (')' * SPECIAL_SQR)
		string = 'fp_exp(' + rewrite_line(string) + ',' + str(2**SPECIAL_SQR) + ')'

		if ' << ' in second:
			# Computing x ^ (2^e)
			e = second.replace(' << ', '')
			if ' + ' not in e:
				e = int(e)
				#string = ('fp_sqr(' * e) + string + (')' * e)
				string = 'fp_exp(' + string + ',' + str(2**e) + ')'
			else:
				assert(' + ' in e)
				if e.count(' + ') == 1:
					e, y = e.split(' + ')
					e = int(e)
					if y == '1': y = '_1'
					#string = 'fp_mul(' + ('fp_sqr(' * e) + string + (')' * e) + ',' + y + ')'
					string = 'fp_mul(fp_exp(' + string + ',' + str(2**e) + ')' + ',' + y + ')'
				else:
					assert(e.count(' + ') == 2)
					e, y, z = e.split(' + ')
					e = int(e)
					if y == '1': y = '_1'
					if z == '1': z = '_1'
					#string = 'fp_mul(fp_mul(' + ('fp_sqr(' * e) + string + (')' * e) + ',' + y + ')' + z + ')'
					string = 'fp_mul(fp_mul(fp_exp(' + string + ',' + str(2**e) + ')' + ',' + y + '),' + z + ')'

		else:
			if (second != ''):
				# Computing x * y
				assert(' << ' not in second)
				assert(' + ' in second)
				if second.count(' + ') == 1:
					y = second.replace(' + ', '')
					if y == '1': y = '_1'
					string = 'fp_mul(' + string + ',' + y + ')'
				else:
					assert(second.count(' + ') == 2)
					xx, y, z = second.split(' + ')
					if y == '1': y = '_1'
					if z == '1': z = '_1'
					string = 'fp_mul(fp_mul(' + string + ',' + y + '),' + z + ')'

		return string

	else:

		# No more required recursion
		assert( ('(' not in line) and (')' not in line) )
		if SPECIAL_SQR:
			assert(' << ' not in line)
			if ' + ' in line:
				if line.count(' + ') == 1:
					x, y = line.split(' + ')
					x = 'fp_sqr(' + x[2:] + ')'
					if y == '1': y = '_1'
					string = 'fp_mul(' + x + ',' + y + ')'
				else:
					assert(line.count(' + ') == 2)
					x, y, z = line.split(' + ')
					x = 'fp_sqr(' + x[2:] + ')'
					if y == '1': y = '_1'
					if z == '1': z = '_1'
					string = 'fp_mul(fp_mul(' + x + ',' + y + '),' + z + ')'
			else:
				string = 'fp_sqr(' + line[2:] + ')'

		elif ' << ' in line:
			# Computing x ^ (2^e)
			x, e = line.split(' << ')
			if x == '1': x = '_1'

			if ' + ' not in e:
				e = int(e)
				#string = ('fp_sqr(' * e) + x + (')' * e)
				string = 'fp_exp(' + x + ',' + str(2**e) + ')'
			else:
				assert(' + ' in e)
				if e.count(' + ') == 1:
					e, y = e.split(' + ')
					if y == '1': y = '_1'
					e = int(e)
					#string = 'fp_mul(' + ('fp_sqr(' * e) + x + (')' * e) + ',' + y + ')'
					string = 'fp_mul(fp_exp(' + x + ',' + str(2**e) + '),' + y + ')'
				else:
					assert(e.count(' + ') == 2)
					e, y, z = e.split(' + ')
					if y == '1': y = '_1'
					if z == '1': z = '_1'
					e = int(e)
					#string = 'fp_mul(fp_mul(' + ('fp_sqr(' * e) + x + (')' * e) + ',' + y + ')' + z + ')'
					string = 'fp_mul(fp_mul(fp_exp(' + x + ',' + str(2**e) + '),' + y + '),'+ z + ')'

		else:
			# Computing x * y
			assert(' << ' not in line)
			assert(' + ' in line)
			if line.count(' + ') == 1:
				x, y = line.split(' + ')
				if x == '1': x = '_1'
				if y == '1': y = '_1'
				string = 'fp_mul(' + x + ',' + y + ')'
			else:
				assert(line.count(' + ') == 2)
				x, y, z = line.split(' + ')
				if x == '1': x = '_1'
				if y == '1': y = '_1'
				if z == '1': z = '_1'
				string = 'fp_mul(fp_mul(' + x + ',' + y + '),' + z + ')'

		return string

def parse_file(filename : str):

	defname = filename.split('/')[-1].replace('.log', '')
	print("from src.fp import fp_mul, fp_sqr, fp_exp")
	print("def %s(x):" % defname)
	print("\t_1 = x")
	dict_addchain = {}
	with open(filename) as FILE:
		for line in FILE:
			if 'return' not in line:
				key, value = line.split(" = ")
			else:
				key = 'return'
				while key in line:
					key = key + ' '

				key = key[:-1]
				value = line.split(key)[1]

			if '2*1' in value:
				value = '2*_1'

			dict_addchain[key] = rewrite_line(value.replace('\n',''))
			print('\t' + key + {True:' ', False:' = '}['return' in key] + dict_addchain[key])

	return dict_addchain

@click.command()
@click.option('--k', default=None, type=str, help="Number of \\ell_i's.")
@click.option('--bits', default=None, type=str, help="Bitlength of the characteristic of the field.")
@click.option('--exponent', default=None, type=click.Choice(("inv", "novem", "quart", "quint", "sept", "sq", "tri")), help="sagemath or magma or raw output format).")
def main(k, bits, exponent):

	original_stdout = sys.stdout # Save a reference to the original standard output
	with open(f'radical/{exponent}_{bits}.py', 'w') as f:
		sys.stdout = f  # Change the standard output to the file we created.
		formula = parse_file(f'crad_primes/p{bits}k{k}_addchains/{exponent}_{bits}.log')

	sys.stdout = original_stdout  # Reset the standard output to its original value
	return formula

if __name__ == "__main__":
	main()