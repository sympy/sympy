# -*- coding: utf-8 -*-

from sympy import *
from matchpy import *
from sympy.integrals.rubi.utility_function import *
from sympy.integrals.rubi.constraints import *
# from sympy.integrals.rubi.symbol import *
from multiset import Multiset
from matchpy.utils import VariableWithCount
from collections import deque
from matchpy.matching.many_to_one import CommutativeMatcher

class CommutativeMatcher2243(CommutativeMatcher):
	_instance = None
	patterns = {
    0: (0, Multiset({}), [
      (VariableWithCount('i3.1.0', 1, 1, None), Mul),
      (VariableWithCount('i3.1.0_1', 1, 1, S(1)), Mul)
])
}
	subjects = {}
	subjects_by_id = {}
	bipartite = BipartiteGraph()
	associative = Mul
	max_optional_count = 1
	anonymous_patterns = set()

	def __init__(self):
		self.add_subject(None)

	@staticmethod
	def get():
		if CommutativeMatcher2243._instance is None:
			CommutativeMatcher2243._instance = CommutativeMatcher2243()
		return CommutativeMatcher2243._instance

	@staticmethod
	def get_match_iter(subject):
		subjects = deque([subject]) if subject is not None else deque()
		subst0 = Substitution()
		# State 2242
		return
		yield


class CommutativeMatcher2239(CommutativeMatcher):
	_instance = None
	patterns = {
    0: (0, Multiset({0: 1}), [
      (VariableWithCount('i3.0', 1, 1, None), Add)
])
}
	subjects = {}
	subjects_by_id = {}
	bipartite = BipartiteGraph()
	associative = Add
	max_optional_count = 0
	anonymous_patterns = set()

	def __init__(self):
		self.add_subject(None)

	@staticmethod
	def get():
		if CommutativeMatcher2239._instance is None:
			CommutativeMatcher2239._instance = CommutativeMatcher2239()
		return CommutativeMatcher2239._instance

	@staticmethod
	def get_match_iter(subject):
		subjects = deque([subject]) if subject is not None else deque()
		subst0 = Substitution()
		# State 2238
		subst1 = Substitution(subst0)
		try:
			subst1.try_add_variable('i3.1.0_1', S(1))
		except ValueError:
			pass
		else:
			# State 2240
			if len(subjects) >= 1:
				tmp2 = subjects.popleft()
				subst2 = Substitution(subst1)
				try:
					subst2.try_add_variable('i3.1.0', tmp2)
				except ValueError:
					pass
				else:
					# State 2241
					if len(subjects) == 0:
						# 0: x*b /; (cons_f2(a, x)) and (cons_f3(b, x)) and (cons_f67(x, a, b))
						yield 0, subst2
				subjects.appendleft(tmp2)
		if len(subjects) >= 1 and isinstance(subjects[0], Mul):
			tmp4 = subjects.popleft()
			associative1 = tmp4
			associative_type1 = type(tmp4)
			subjects5 = deque(tmp4._args)
			matcher = CommutativeMatcher2243.get()
			tmp6 = subjects5
			subjects5 = []
			for s in tmp6:
				matcher.add_subject(s)
			for pattern_index, subst1 in matcher.match(tmp6, subst0):
				if pattern_index == 0:
					# State 2244
					if len(subjects) == 0:
						# 0: x*b /; (cons_f2(a, x)) and (cons_f3(b, x)) and (cons_f67(x, a, b))
						yield 0, subst1
			subjects.appendleft(tmp4)
		return
		yield

def match_root(subject):
	subjects = deque([subject]) if subject is not None else deque()
	subst0 = Substitution()
	# State 2222
	if len(subjects) >= 1 and isinstance(subjects[0], Integral):
		tmp1 = subjects.popleft()
		subjects2 = deque((tmp1._args[0],) + tmp1._args[1])
		# State 2223
		if len(subjects2) >= 1 and isinstance(subjects2[0], Pow):
			tmp3 = subjects2.popleft()
			subjects4 = deque(tmp3._args)
			# State 2224
			if len(subjects4) >= 1:
				tmp5 = subjects4.popleft()
				subst1 = Substitution(subst0)
				try:
					subst1.try_add_variable('i2', tmp5)
				except ValueError:
					pass
				else:
					if not ('i2' in subst1 and 'i3' in subst1) or cons_f21(subst1['i3'], subst1['i2']):
						# State 2225
						if len(subjects4) >= 1 and subjects4[0] == -1:
							tmp7 = subjects4.popleft()
							# State 2226
							if len(subjects4) == 0:
								# State 2227
								if len(subjects2) >= 1:
									tmp8 = subjects2.popleft()
									subst2 = Substitution(subst1)
									try:
										subst2.try_add_variable('i2', tmp8)
									except ValueError:
										pass
									else:
										# State 2228
										if len(subjects2) == 0:
											# State 2229
											if len(subjects) == 0:
												pass
									subjects2.appendleft(tmp8)
							subjects4.appendleft(tmp7)
						subst2 = Substitution(subst1)
						try:
							subst2.try_add_variable('i3', 1)
						except ValueError:
							pass
						else:
							if not ('i2' in subst2 and 'i3' in subst2) or cons_f21(subst2['i3'], subst2['i2']):
								if not ('i3' in subst2) or cons_f66(subst2['i3']):
									# State 2234
									if len(subjects4) == 0:
										# State 2235
										if len(subjects2) >= 1:
											tmp11 = subjects2.popleft()
											subst3 = Substitution(subst2)
											try:
												subst3.try_add_variable('i2', tmp11)
											except ValueError:
												pass
											else:
												if not ('i2' in subst3 and 'i3' in subst3) or cons_f21(subst3['i3'], subst3['i2']):
													# State 2236
													if len(subjects2) == 0:
														# State 2237
														if len(subjects) == 0:
															tmp_subst = Substitution()
															tmp_subst['x'] = subst3['i2']
															tmp_subst['m'] = subst3['i3']
															# 1: Integral(x**m, x) /; (cons_f21(m, x)) and (cons_f66(m))
															yield 38, tmp_subst
											subjects2.appendleft(tmp11)
						if len(subjects4) >= 1:
							tmp13 = subjects4.popleft()
							subst2 = Substitution(subst1)
							try:
								subst2.try_add_variable('i3', tmp13)
							except ValueError:
								pass
							else:
								if not ('i2' in subst2 and 'i3' in subst2) or cons_f21(subst2['i3'], subst2['i2']):
									if not ('i3' in subst2) or cons_f66(subst2['i3']):
										# State 2234
										if len(subjects4) == 0:
											# State 2235
											if len(subjects2) >= 1:
												tmp15 = subjects2.popleft()
												subst3 = Substitution(subst2)
												try:
													subst3.try_add_variable('i2', tmp15)
												except ValueError:
													pass
												else:
													if not ('i2' in subst3 and 'i3' in subst3) or cons_f21(subst3['i3'], subst3['i2']):
														# State 2236
														if len(subjects2) == 0:
															# State 2237
															if len(subjects) == 0:
																tmp_subst = Substitution()
																tmp_subst['x'] = subst3['i2']
																tmp_subst['m'] = subst3['i3']
																# 1: Integral(x**m, x) /; (cons_f21(m, x)) and (cons_f66(m))
																yield 38, tmp_subst
												subjects2.appendleft(tmp15)
							subjects4.appendleft(tmp13)
					# State 2225
					if len(subjects4) >= 1 and subjects4[0] == -1:
						tmp17 = subjects4.popleft()
						# State 2226
						if len(subjects4) == 0:
							# State 2227
							if len(subjects2) >= 1:
								tmp18 = subjects2.popleft()
								subst2 = Substitution(subst1)
								try:
									subst2.try_add_variable('i2', tmp18)
								except ValueError:
									pass
								else:
									# State 2228
									if len(subjects2) == 0:
										# State 2229
										if len(subjects) == 0:
											tmp_subst = Substitution()
											tmp_subst['x'] = subst2['i2']
											# 0: Integral(1/x, x)
											yield 37, tmp_subst
								subjects2.appendleft(tmp18)
						subjects4.appendleft(tmp17)
					subst2 = Substitution(subst1)
					try:
						subst2.try_add_variable('i3', 1)
					except ValueError:
						pass
					else:
						pass
					if len(subjects4) >= 1:
						tmp21 = subjects4.popleft()
						subst2 = Substitution(subst1)
						try:
							subst2.try_add_variable('i3', tmp21)
						except ValueError:
							pass
						else:
							pass
						subjects4.appendleft(tmp21)
				subjects4.appendleft(tmp5)
			if len(subjects4) >= 1 and isinstance(subjects4[0], Add):
				tmp23 = subjects4.popleft()
				associative1 = tmp23
				associative_type1 = type(tmp23)
				subjects24 = deque(tmp23._args)
				matcher = CommutativeMatcher2239.get()
				tmp25 = subjects24
				subjects24 = []
				for s in tmp25:
					matcher.add_subject(s)
				for pattern_index, subst1 in matcher.match(tmp25, subst0):
					if pattern_index == 0:
						pass
				subjects4.appendleft(tmp23)
			subjects2.appendleft(tmp3)
		subst1 = Substitution(subst0)
		try:
			subst1.try_add_variable('i3', S(1))
		except ValueError:
			pass
		else:
			# State 2230
			if len(subjects2) >= 1:
				tmp27 = subjects2.popleft()
				subst2 = Substitution(subst1)
				try:
					subst2.try_add_variable('i2', tmp27)
				except ValueError:
					pass
				else:
					if not ('i2' in subst2 and 'i3' in subst2) or cons_f21(subst2['i3'], subst2['i2']):
						# State 2231
						if len(subjects2) >= 1:
							tmp29 = subjects2.popleft()
							subst3 = Substitution(subst2)
							try:
								subst3.try_add_variable('i2', tmp29)
							except ValueError:
								pass
							else:
								if not ('i2' in subst3 and 'i3' in subst3) or cons_f21(subst3['i3'], subst3['i2']):
									# State 2232
									if len(subjects2) == 0:
										# State 2233
										if len(subjects) == 0:
											tmp_subst = Substitution()
											tmp_subst['x'] = subst3['i2']
											tmp_subst['m'] = subst3['i3']
											# 1: Integral(x**m, x) /; (cons_f21(m, x)) and (cons_f66(m))
											yield 38, tmp_subst
							subjects2.appendleft(tmp29)
				subjects2.appendleft(tmp27)
		subjects.appendleft(tmp1)
	return
	yield
