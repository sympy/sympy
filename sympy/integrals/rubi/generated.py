# -*- coding: utf-8 -*-

from sympy import *
from matchpy import *
from sympy.integrals.rubi.utility_function import *
from sympy.integrals.rubi.constraints import *
# from sympy.integrals.rubi.symbol import *
from collections import deque
def match_root(subject):
	subjects = deque([subject]) if subject is not None else deque()
	subst0 = Substitution()
	# State 2210
	if len(subjects) >= 1 and isinstance(subjects[0], Integral):
		tmp1 = subjects.popleft()
		subjects2 = deque((tmp1._args[0],) + tmp1._args[1])
		# State 2211
		if len(subjects2) >= 1 and isinstance(subjects2[0], Pow):
			tmp3 = subjects2.popleft()
			subjects4 = deque(tmp3._args)
			# State 2212
			if len(subjects4) >= 1:
				tmp5 = subjects4.popleft()
				subst1 = Substitution(subst0)
				try:
					subst1.try_add_variable('i2', tmp5)
				except ValueError:
					pass
				else:
					if 'i3' in subst1 and 'i2' in subst1 and cons_f21(subst1['i3'], subst1['i2']):
						# State 2213
						if len(subjects4) >= 1 and subjects4[0] == -1:
							tmp7 = subjects4.popleft()
							# State 2214
							if len(subjects4) == 0:
								# State 2215
								if len(subjects2) >= 1:
									tmp8 = subjects2.popleft()
									subst2 = Substitution(subst1)
									try:
										subst2.try_add_variable('i2', tmp8)
									except ValueError:
										pass
									else:
										# State 2216
										if len(subjects2) == 0:
											# State 2217
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
							if 'i3' in subst2 and 'i2' in subst2 and cons_f21(subst2['i3'], subst2['i2']):
								if 'i3' in subst2 and cons_f66(subst2['i3']):
									# State 2222
									if len(subjects4) == 0:
										# State 2223
										if len(subjects2) >= 1:
											tmp11 = subjects2.popleft()
											subst3 = Substitution(subst2)
											try:
												subst3.try_add_variable('i2', tmp11)
											except ValueError:
												pass
											else:
												if 'i3' in subst3 and 'i2' in subst3 and cons_f21(subst3['i3'], subst3['i2']):
													# State 2224
													if len(subjects2) == 0:
														# State 2225
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
								if 'i3' in subst2 and 'i2' in subst2 and cons_f21(subst2['i3'], subst2['i2']):
									if 'i3' in subst2 and cons_f66(subst2['i3']):
										# State 2222
										if len(subjects4) == 0:
											# State 2223
											if len(subjects2) >= 1:
												tmp15 = subjects2.popleft()
												subst3 = Substitution(subst2)
												try:
													subst3.try_add_variable('i2', tmp15)
												except ValueError:
													pass
												else:
													if 'i3' in subst3 and 'i2' in subst3 and cons_f21(subst3['i3'], subst3['i2']):
														# State 2224
														if len(subjects2) == 0:
															# State 2225
															if len(subjects) == 0:
																tmp_subst = Substitution()
																tmp_subst['x'] = subst3['i2']
																tmp_subst['m'] = subst3['i3']
																# 1: Integral(x**m, x) /; (cons_f21(m, x)) and (cons_f66(m))
																yield 38, tmp_subst
												subjects2.appendleft(tmp15)
							subjects4.appendleft(tmp13)
					# State 2213
					if len(subjects4) >= 1 and subjects4[0] == -1:
						tmp17 = subjects4.popleft()
						# State 2214
						if len(subjects4) == 0:
							# State 2215
							if len(subjects2) >= 1:
								tmp18 = subjects2.popleft()
								subst2 = Substitution(subst1)
								try:
									subst2.try_add_variable('i2', tmp18)
								except ValueError:
									pass
								else:
									# State 2216
									if len(subjects2) == 0:
										# State 2217
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
			subjects2.appendleft(tmp3)
		subst1 = Substitution(subst0)
		try:
			subst1.try_add_variable('i3', S(1))
		except ValueError:
			pass
		else:
			# State 2218
			if len(subjects2) >= 1:
				tmp24 = subjects2.popleft()
				subst2 = Substitution(subst1)
				try:
					subst2.try_add_variable('i2', tmp24)
				except ValueError:
					pass
				else:
					if 'i3' in subst2 and 'i2' in subst2 and cons_f21(subst2['i3'], subst2['i2']):
						# State 2219
						if len(subjects2) >= 1:
							tmp26 = subjects2.popleft()
							subst3 = Substitution(subst2)
							try:
								subst3.try_add_variable('i2', tmp26)
							except ValueError:
								pass
							else:
								if 'i3' in subst3 and 'i2' in subst3 and cons_f21(subst3['i3'], subst3['i2']):
									# State 2220
									if len(subjects2) == 0:
										# State 2221
										if len(subjects) == 0:
											tmp_subst = Substitution()
											tmp_subst['x'] = subst3['i2']
											tmp_subst['m'] = subst3['i3']
											# 1: Integral(x**m, x) /; (cons_f21(m, x)) and (cons_f66(m))
											yield 38, tmp_subst
							subjects2.appendleft(tmp26)
				subjects2.appendleft(tmp24)
		subjects.appendleft(tmp1)
	return
	yield
