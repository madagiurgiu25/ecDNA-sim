import os
import pprint
import random

random.seed(12345)
CHRLEN = {"chr2": 242193529,
          "chr8": 145138636,
          "chr12": 133275309}

# allowed chromosomal regions
CHRRANGES = {"chr2": (4880899, 83143574),
             "chr8": (50000000, 106259331),
             "chr12": (48581298, 120529644)}

# n = number of fragments
# single_chr = intrachromosomal (True), interchromosomal (False)

N = "n"
NEIGHBOR = "neighbor"
MIN_TANDEM_DUPLICATION = "min_tandem_duplication"
SINGLE_CHR = "single_chr"
WEIGHT_CHR = "weights_chr"
P_DEL_LEFT = "p_del_left"
P_DEL_RIGHT = "p_del_right"
P_INVERT = "p_invert"
P_DUP = "p_dup"
P_FOLDBACK = "p_foldback"
P_RETURN = "p_return"
P_CHOOSE_RANDOM = "p_choose_random"

# variables included in the metainformation string
N_FRAG = "fragments"
FRAG_LEN = "frag_len"
READ_LEN = "read_len"
SMALL_DEL = "small_del"
DUP = "dup"
INV = "inv"
INTERCHR = "interchromosomal"
MULTI_REGION = "multi-region"
FOLDBACK = "foldback"
RETURN = "return"  # new fragment is next to the other

N_FRAG_POS = 0
FRAG_LEN_POS = 1
# READ_LEN_POS = 1
SMALL_DEL_POS = 2
DUP_POS = 3
INV_POS = 4
INTERCHR_POS = 5
MULTI_REGION_POS = 6
FOLDBACK_POS = 7
RETURN_POS = 8

CHR = 0
START = 1
END = 2
STRAND = 3

EVENTS = {
	N_FRAG: N_FRAG_POS,
	FRAG_LEN: FRAG_LEN_POS,
	SMALL_DEL: SMALL_DEL_POS,
	DUP: DUP_POS,
	INV: INV_POS,
	INTERCHR: INTERCHR_POS,
	MULTI_REGION: MULTI_REGION_POS,
	FOLDBACK: FOLDBACK_POS,
	RETURN: RETURN_POS
}

CONF_INIT = {
	N: 3,
	NEIGHBOR: False,
	SINGLE_CHR: True,
	MIN_TANDEM_DUPLICATION: 1000,
	WEIGHT_CHR: [50, 50],
	P_DEL_LEFT: 20,
	P_DEL_RIGHT: 30,
	P_INVERT: 40,
	P_DUP: 10,
	P_FOLDBACK: 5,
	P_RETURN: 5,
	# probl to choose random position or otherwise transit to foldback or return
	P_CHOOSE_RANDOM: 30
}

pp = pprint.PrettyPrinter(indent=4)


def generate_binary(length, l, all):
	if length > 0:
		generate_binary(length - 1, l + [0], all)
		generate_binary(length - 1, l + [1], all)
	else:
		all.append(l)


def generate_conformation():
	"""
	Generate all possible conformation which we want to simulate
	"""
	# N_FRAG_POS = 0
	# FRAG_LEN_POS = 1
	# SMALL_DEL_POS = 2
	# DUP_POS = 3
	# INV_POS = 4
	# INTERCHR_POS = 5
	# MULTI_REGION_POS = 6
	# FOLDBACK_POS = 7
	# RETURN_POS = 8

	# initialize simple excision
	all_conf = [[1, 0, 0, 0, 0, 0, 0, 0, 0]]

	values = {N_FRAG: list(range(1, 10)),
	          SMALL_DEL: [0, 1],
	          DUP: [0, 1],
	          INV: [0, 1],
	          INTERCHR: [0, 1],
	          MULTI_REGION: [0, 1],
	          FOLDBACK: [0, 1],
	          RETURN: [0, 1]
	          }

	bins = []
	generate_binary(6, [], bins)

	for f in range(2, 11):
		for b in bins:
			# [#frag, len] + [small_del, dup, inv, interch, multi+region, foldback] + [return]
			newconf = [f,0]+b+[0]
			all_conf.append(newconf)

	return all_conf


def code_metainformation(topology):
	"""
	Configuration string for every ecDNA
	Convert configuration dictionary to a meta information array.

	Args:
		topology (dict): Contains all the configuration values
	"""
	# array of information
	# pos 0 - read_len {0,1,2,3} - 0 - << 7000bp, 1 - ~7000, 2 - >> 7000, 3
	#
	arr = [0] * 9
	arr[DUP_POS] = 1 if topology[DUP] > 0 else 0
	arr[INV_POS] = 1 if topology[INV] > 0 else 0
	arr[INTERCHR_POS] = 1 if topology[INTERCHR] > 0 else 0
	arr[MULTI_REGION_POS] = 1 if topology[MULTI_REGION] > 0 else 0
	arr[FOLDBACK_POS] = 1 if topology[FOLDBACK] > 0 else 0
	return "".join([str(v) for v in arr])


def is_multiregion(e1, e2, threshold=5000):
	"""
	Consider multiregion if the distance between fragments is min >= threshold
	"""
	# sit on different chromosomes
	if e1[CHR] != e2[CHR]:
		return True
	if e1[END] < e2[START] and e2[START] - e1[END] >= threshold:
		return True
	if e1[START] > e2[END] and e1[START] - e2[END] >= threshold:
		return True
	return False


def pick_chr(chrlen, conf):
	"""
	Choose ecDNA origin chromosome
	"""
	return random.choices([p for p in chrlen], weights=conf[WEIGHT_CHR], k=1)[0]


def pick_pos(chrranges, conf, chr_):
	"""
	Pick genomic position on chr
	"""
	return random.randrange(chrranges[chr_][0], chrranges[chr_][1], 5)


def pick_pos_overlap(ecDNA, index):
	"""
	Pick genomic position that overlaps the previous fragment (start,stop).

	Args:
		ecDNA (list): List of fragments
		index (int): Index of the previous neighboring element
	"""
	start = min(ecDNA[index][START], ecDNA[index][END])
	end = max(ecDNA[index][START], ecDNA[index][END])
	return random.randrange(start, end, 5)


def pick_len(chrlen, conf, chr, start, meansize=300000, sd=5000):
	"""
	Pick length of ecDNA fragment
	"""
	return random.randrange(meansize - sd, meansize + sd, 100)


def pick_len_overlap(ecDNA, previous_index):
	"""
	Pick length of ecDNA fragment
	"""
	previous_len = abs(ecDNA[previous_index][START] - ecDNA[previous_index][END])
	return random.randrange(int(previous_len / 4), int(previous_len / 2), 100)


def choose_fragment_random(chrlen, chrranges, conf):
	"""
	Choose genomic fragment of origin
	"""
	chr_ = pick_chr(chrlen, conf)
	start_ = pick_pos(chrranges, conf, chr_)
	len_ = pick_len(chrlen, chrranges, chr_, start_)
	end_ = start_ + len_

	return chr_, start_, end_, len_


def choose_fragment_neighbor(chrlen, chrranges, conf, index, ecDNA):
	"""
	Choose genomic fragment of origin
	"""
	if index == -1:
		chr_ = pick_chr(chrlen, conf)
		start_ = pick_pos(chrranges, conf, chr_)
		len_ = pick_len(chrlen, conf, chr_, start_)
		end_ = start_ + len_
	else:
		# 0 1 2 3 4 5
		# chr_,start_,end_,strand_,count_frag
		chr_ = ecDNA[index][0]
		start_ = ecDNA[index][2]
		len_ = pick_len(chrlen, conf, chr_, start_)
		end_ = start_ + len_

	return chr_, start_, end_, len_


def choose_fragment_pseudorandom(chrlen, chrranges, conf, index, ecDNA):
	"""
	Choose pseudorandom only with a probability is random
	"""
	neighbor = random.choices([1, 0], k=1)[0]

	# choose to pick a neighboring fragment
	if neighbor == 1:
		chr_, start_, end_, len_ = choose_fragment_neighbor(chrlen, chrranges, conf, index, ecDNA)
	else:
		# choose random position on the same chromosome
		if conf[SINGLE_CHR]:
			# initialize chromosome
			chr_ = pick_chr(chrlen, conf) if index == -1 else ecDNA[index][0]
		else:
			# pick random chromosome
			chr_ = pick_chr(chrlen, conf)

		start_ = pick_pos(chrranges, conf, chr_)
		len_ = pick_len(chrlen, conf, chr_, start_)
		end_ = start_ + len_

	return chr_, start_, end_, len_


def choose_fragment(chrlen, chrranges, index, conf, ecDNA, topology, type="random"):
	if conf[NEIGHBOR]:
		chr_, start_, end_, len_ = choose_fragment_neighbor(chrlen, chrranges, conf, index, ecDNA)
	elif type == "pseudorandom":
		chr_, start_, end_, len_ = choose_fragment_pseudorandom(chrlen, chrranges, conf, index, ecDNA)
	else:
		chr_, start_, end_, len_ = choose_fragment_random(chrlen, chrranges, conf)

	return chr_, start_, end_, len_


def choose_fold(chrlen, chrranges, previous_index, conf, ecDNA, topology):
	"""
	Choose foldback fragment
	"""
	chr_ = ecDNA[previous_index][0]
	start_ = pick_pos_overlap(ecDNA, previous_index)
	len_ = pick_len_overlap(ecDNA, previous_index)
	end_ = start_ + len_

	return chr_, start_, end_, len_


def crop_left_end(fragment_len, p):
	"""
	Decide to crop (10% of fragment len) or not the left end
	"""
	del_left = random.choices([1, 0], weights=[p, 100 - p], k=1)[0]
	return del_left, 0.1 * fragment_len


def crop_right_end(fragment_len, p):
	"""
	Decide to crop (10% of fragment len) or not the right end
	"""
	del_right = random.choices([1, 0], weights=[p, 100 - p], k=1)[0]
	return del_right, 0.1 * fragment_len


def small_deletions(chrlen, chr_, start_, end_, len_, conf, topology):
	"""
	Emulate small deletions
	"""

	(left, lsize) = crop_left_end(len_, conf[P_DEL_LEFT])
	(right, rsize) = crop_right_end(len_, conf[P_DEL_RIGHT])

	# crop left part
	if left == 1:
		start_ = min(chrlen[chr_], start_ + lsize)
		topology[SMALL_DEL] += 1

	# crop right part
	if right == 1:
		end_ = max(0, end_ - rsize)
		topology[SMALL_DEL] += 1

	return int(start_), int(end_)


def invert_event(conf):
	"""
	Invert fragment with a probability conf[P_INVERT]
	"""
	p = conf[P_INVERT]
	return random.choices([1, 0], weights=[p, 100 - p], k=1)[0]


def duplicate_event(conf):
	"""
	Duplicate fragment with a probability conf[P_DUP]
	"""
	p = conf[P_DUP]
	return random.choices([1, 0], weights=[p, 100 - p], k=1)[0]


def invert(conf, topology):
	strand = "+"
	if invert_event(conf) == 1:
		topology[INV] += 1
		strand = "-"
	return strand


def duplicate(conf, topology):
	if duplicate_event(conf) == 1:
		topology[DUP] += 1
		return True
	return False


def pick_nr_copies(conf):
	"""
	Pick number of ecDNA copies
	"""
	return random.choices(range(400, 700, 10), k=1)[0]


def foldback(conf, topology):
	"""
	Foldback fragment with a probability conf[P_FOLDBACK]
	"""
	p = conf[P_FOLDBACK]
	return random.choices([1, 0], weights=[p, 100 - p], k=1)[0]


def back(conf, topology):
	"""
	Backstructure conf[P_RETURN]
	"""
	p = conf[P_FOLDBACK]
	return random.choices([1, 0], weights=[p, 100 - p], k=1)[0]


def choose(conf):
	"""
	0 - Random position
	1 - Foldback
	2 - Return
	"""
	p = conf[P_CHOOSE_RANDOM]
	return random.choices([0, 1, 2], weights=[p, 3 * (100 - p) / 4, (100 - p) / 4], k=1)[0]


def compare_metainformation(acceptance_criteria, metastring):
	for key in metastring:
		if key in EVENTS:
			id = EVENTS[key]
			if acceptance_criteria[id] == 1 and metastring[key] == 0:
				return False
	return True


def simulation(chrlen, chrranges, conf, type="random", acceptance_criteria=[3, 0, 1, 1, 0, 1, 1, 1, 0]):
	"""
	Simulate a single ecDNA structure
	"""

	# 1. add all default parameters (unchanged) to conf file
	for key in CONF_INIT:
		if key not in conf:
			conf[key] = CONF_INIT[key]

	# pp.pprint(conf)

	pass_structure = False
	max_iterations = 300
	count_iter = 0

	while (not pass_structure) and count_iter < max_iterations:

		ecDNA = []
		count_frag = -1
		count_iter += 1

		topology = {N_FRAG: 0,
		            FRAG_LEN: 0,
		            SMALL_DEL: 0,
		            DUP: 0,
		            INV: 0,
		            INTERCHR: 0,
		            MULTI_REGION: 0,
		            FOLDBACK: 0,
		            RETURN: 0
		            }

		# 2. start simulation
		while len(ecDNA) < conf[N]:

			eventtype = choose(conf)

			if eventtype != 0 and count_frag == -1:
				continue

			if eventtype == 0:
				# 2.1. choose fragment
				chr_, start_, end_, len_ = choose_fragment(chrlen, chrranges, count_frag, conf, ecDNA, topology,
				                                           type=type)
			# print("Choose fragment Nr. chr start end len, ", count_frag, chr_, start_, end_, len_)
			elif eventtype == 1:
				# foldback
				topology[FOLDBACK] += 1
				chr_, start_, end_, len_ = choose_fold(chrlen, chrranges, count_frag, conf, ecDNA, topology)
			else:
				# return/back
				pass

			# 2.2 emulate small deletions
			start_, end_ = small_deletions(chrlen, chr_, start_, end_, len_, conf, topology)
			# print("Small del Nr. chr start end len, ", count_frag, chr_, start_, end_, len_)

			# 2.3 emulate inversion
			strand_ = invert(conf, topology)

			count_frag += 1
			topology[N_FRAG] += 1
			ecDNA.append((chr_, start_, end_, strand_, count_frag))

			# mark topology as interchromosomal
			if ecDNA[count_frag - 1][0] != ecDNA[count_frag][0]:
				topology[INTERCHR] += 1

			# multiregion
			if is_multiregion(ecDNA[count_frag - 1], ecDNA[count_frag]):
				topology[MULTI_REGION] += 1

			# 2.4 allow simple duplication (1>kbp)
			if abs(start_ - end_) > conf[MIN_TANDEM_DUPLICATION]:
				while ((count_frag + 1) < conf[N]) and duplicate(conf, topology):
					# if count_frag < conf[N]:
					#     dup_ = duplicate(conf)
					#     if dup_ == 1:
					count_frag += 1
					topology[N_FRAG] += 1
					ecDNA.append((chr_, start_, end_, strand_, count_frag))

		# 3. choose coverage
		cov_ = pick_nr_copies(conf)

		# 4. check if conformation correct

		pass_structure = compare_metainformation(acceptance_criteria, topology)

	# if pass_structure:
	# 	print("Pass")
	# 	print(acceptance_criteria)
	# 	print(topology)
	# 	print(ecDNA)
	# 	print()
	# else:
	# 	print("No pass")
	# 	print(ecDNA)
	# 	print(topology)
	# 	print()

	if count_iter == max_iterations:
		return None, None, topology

	return ecDNA, cov_, topology


def ecDNA2bed(ecDNAdict, fout):
	# print(os.getcwd())
	# print(os.listdir())

	with open(fout, "w") as f:
		f.write("#chr\tstart\tstop\tdirection\ttarget\tcoverage\tstructure\tfragment\n")
		for circ in ecDNAdict:
			for i, frag in enumerate(ecDNAdict[circ]["structure"]):
				chr_, start_, end_, strand_, count_frag = frag
				strand_ = "+" if strand_ == 0 else "-"
				coverage_ = ecDNAdict[circ]["coverage"]
				circ_ = """circ_{}""".format(str(circ))
				circ_ = """circ_{}""".format(str(circ))
				f.write("""{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n""".format(chr_,
				                                                      str(start_),
				                                                      str(end_),
				                                                      strand_,
				                                                      circ_,
				                                                      str(coverage_),
				                                                      circ_,
				                                                      str(i)
				                                                      ))

def ecDNA2bed_single(ecDNAdict, fout):
	# print(os.getcwd())
	# print(os.listdir())

	os.makedirs(os.path.dirname(fout), exist_ok=True)
	with open(fout, "w") as f:
		f.write("#chr\tstart\tstop\tdirection\ttarget\tcoverage\tstructure\tfragment\n")

		for frag in ecDNAdict["structure"]:
			chr_, start_, end_, strand_, count_frag = frag
			strand_ = strand_
			coverage_ = ecDNAdict["coverage"]
			circ_ = """circ_0"""
			f.write("""{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n""".format(chr_,
			                                                      str(start_),
			                                                      str(end_),
			                                                      strand_,
			                                                      circ_,
			                                                      str(coverage_),
			                                                      circ_,
			                                                      str(count_frag)
			                                                      ))
