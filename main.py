import json
import pprint
from collections import defaultdict

import simulate as s

pp = pprint.PrettyPrinter(indent=4)


def simulate_simple_excisions(chrlen, chrranges):
	"""
	Simulate simple excision.
	"""
	conf = {s.N: 1,
	        s.NEIGHBOR: True,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 0,
	        s.P_DEL_RIGHT: 0,
	        s.P_INVERT: 0,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        s.P_CHOOSE_RANDOM: 1
	        }

	for i in range(0, 5):
		ecDNAdict = defaultdict(int)
		file = """sim_2/excision/circ_1_{}.bed""".format(str(i))
		# print(file)
		for j in range(0, conf[s.N]):
			ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="random")
			ecDNAdict[j] = {}
			ecDNAdict[j]["structure"] = ecDNA
			ecDNAdict[j]["coverage"] = cov
		s.ecDNA2bed(ecDNAdict, file)


def simulate_simple_deletions(chrlen, chrranges):
	"""
	Simulate simple events like small deletions
	"""
	conf = {s.N: 2,
	        s.NEIGHBOR: True,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 0,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="nonrandom")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/simple_del/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_simple_inversion(chrlen, chrranges):
	"""
	Simulate simple inversion and deletions
	"""
	conf = {s.N: 2,
	        s.NEIGHBOR: True,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 0,
	        s.P_DEL_RIGHT: 0,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="nonrandom")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/simple_inv/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_intrachrom_multi_region_random(chrlen, chrranges):
	"""
	Simulate intrachrom multi region
	"""
	conf = {s.NEIGHBOR: False,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        # related to foldback
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="random")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/intrachrom_random/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_intrachrom_multi_region_pseudorandom(chrlen, chrranges):
	"""
	Simulate intrachrom multi region
	"""
	conf = {s.NEIGHBOR: False,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        # related to foldbacks
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="pseudorandom")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/intrachrom/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_simple_mix(chrlen, chrranges):
	"""
	Simulate simple inversion and deletions
	"""
	conf = {s.NEIGHBOR: True,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        # related to foldbacks
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="random")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/simple_mix/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_interchr_multi_regions(chrlen, chrranges):
	conf = {s.NEIGHBOR: False,
	        s.SINGLE_CHR: False,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 0,
	        s.P_RETURN: 0,
	        # related to foldbacks
	        s.P_CHOOSE_RANDOM: 1
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf, type="pseudorandom")
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/interchrom/circ_{}_{}.bed""".format(str(k), str(i)))


def simulate_foldbacks(chrlen, chrranges):
	conf = {s.NEIGHBOR: True,
	        s.SINGLE_CHR: True,
	        s.MIN_TANDEM_DUPLICATION: 1000,
	        s.WEIGHT_CHR: [50] * len(chrlen),
	        s.P_DEL_LEFT: 60,
	        s.P_DEL_RIGHT: 60,
	        s.P_INVERT: 40,
	        s.P_DUP: 0,
	        s.P_FOLDBACK: 40,
	        s.P_RETURN: 0,
	        # related to foldbacks
	        s.P_CHOOSE_RANDOM: 10
	        }

	# define number of fragments
	for k in range(2, 10):
		conf[s.N] = k

		# simulation number
		for i in [1, 2, 3]:
			ecDNAdict = defaultdict(int)

			# how many ecDNA per file
			for j in range(0, 1):
				ecDNA, cov = s.simulation(chrlen, chrranges, conf)
				ecDNAdict[j] = {}
				ecDNAdict[j]["structure"] = ecDNA
				ecDNAdict[j]["coverage"] = cov
			s.ecDNA2bed(ecDNAdict, """sim_2/foldback/circ_{}_{}.bed""".format(str(k), str(i)))


def generate_allowed(specify_conformation=None):
	if not specify_conformation:
		return s.generate_conformation()
	else:
		return [specify_conformation]


def simulate_structure(chrlen, chrranges, repeat_simulation=5, specify_conformation=[3, 0, 1, 1, 0, 1, 1, 1, 0]):
	"""
	Args:
		chrlen (dict):
		chrranges (dict):
		accept_criteria (list): topology = {N_FRAG: 0,
	            FRAG_LEN: 0,
	            SMALL_DEL: 0,
	            DUP: 0,
	            INV: 0,
	            INTERCHR: 0,
	            MULTI_REGION: 0,
	            FOLDBACK: 0,
	            RETURN: 0
	            }
	"""
	allowed_conformations = generate_allowed(specify_conformation)
	count = 0
	batch = 0
	with open("sim_all/summary.txt", "w") as f:

		f.write("#sim_id\tbatch_id\tcount_frag\tsmall_del\tdup\tinv\tinterchr\tmulti_region\tfoldback\ttopology\n")

		# simulate
		for j, accept_criteria in enumerate(allowed_conformations):
			conf = {}
			conf[s.N] = accept_criteria[s.N_FRAG_POS]
			conf[s.FRAG_LEN] = accept_criteria[s.FRAG_LEN_POS]
			conf[s.NEIGHBOR] = True if accept_criteria[s.MULTI_REGION_POS] == 0 and accept_criteria[
				s.INTERCHR_POS] == 0 else False
			conf[s.MIN_TANDEM_DUPLICATION] = 1000
			conf[s.WEIGHT_CHR] = [50] * len(chrlen)
			conf[s.P_DEL_LEFT] = 0 if accept_criteria[s.SMALL_DEL_POS] == 0 else 30
			conf[s.P_DEL_RIGHT] = 0 if accept_criteria[s.SMALL_DEL_POS] == 0 else 30
			conf[s.P_DUP] = 0 if accept_criteria[s.DUP_POS] == 0 else 30
			conf[s.P_INVERT] = 0 if accept_criteria[s.INV_POS] == 0 else 30
			conf[s.P_FOLDBACK] = 0 if accept_criteria[s.FOLDBACK_POS] == 0 else 30
			conf[s.P_RETURN] = 0
			conf[s.P_CHOOSE_RANDOM] = 100 if accept_criteria[s.FOLDBACK_POS] == 0 and accept_criteria[
				s.RETURN_POS] == 0 else 20
			conf["chrs_len"] = chrlen
			conf["chrs_ranges"] = chrranges
			conf["accept_criteria_description"] = "[#frag, len, small_del, dup, inv, interch, multi+region, foldback, return]"

			for i in range(0, repeat_simulation):

				conf["sim_iteration"] = i
				conf["accept_criteria"] = accept_criteria

				ecDNA, cov, topology = s.simulation(chrlen, chrranges, conf, acceptance_criteria=accept_criteria)

				if ecDNA == None:
					# print("Warning! could not find topology", accept_criteria)
					pass
				else:
					ecDNAdict = {}
					ecDNAdict["structure"] = ecDNA
					ecDNAdict["coverage"] = cov
					ecDNAdict["topology"] = topology
					metainformation = s.code_metainformation(topology)
					count += 1
					prefixfile = """sim_all/batch_{}/sim_{}""".format(str(batch), str(count))
					f.write("""sim_{}\tbatch_{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n""".format(str(count),
					                                                                    str(batch),
					                                                                    topology[s.N_FRAG],
					                                                                    topology[s.DUP],
					                                                                    topology[s.INV],
					                                                                    topology[s.INTERCHR],
					                                                                    topology[s.MULTI_REGION],
					                                                                    topology[s.FOLDBACK],
					                                                                    metainformation))
					s.ecDNA2bed_single(ecDNAdict, prefixfile + ".bed")
					with open(prefixfile + ".json", 'w') as g:
						json.dump(conf, g)
					if count % 100 == 0:
						batch += 1

	print("#Simulation templates ", count)
