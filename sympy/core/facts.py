# -*- coding: utf-8 -*-

"""This is rule-based deduction system for SymPy

The whole thing is split into two parts

 - rules compilation and preparation of tables
 - runtime inference

For rule-based inference engines, the classical work is RETE algorithm [1],
[2] Although we are not implementing it in full (or even significantly)
it's still still worth a read to understand the underlying ideas.

In short, every rule in a system of rules is one of two forms:

 - atom                     -> ...      (alpha rule)
 - And(atom1, atom2, ...)   -> ...      (beta rule)


The major complexity is in efficient beta-rules processing and usually for an
expert system a lot of effort goes into code that operates on beta-rules.


Here we take minimalistic approach to get something usable first.

 - (preparation)    of alpha- and beta- networks, everything except
 - (runtime)        FactRules.deduce_all_facts

             _____________________________________
            ( Kirr: I've never thought that doing )
            ( logic stuff is that difficult...    )
             -------------------------------------
                    o   ^__^
                     o  (oo)\_______
                        (__)\       )\/\
                            ||----w |
                            ||     ||


Some references on the topic
----------------------------

[1] http://en.wikipedia.org/wiki/Rete_algorithm
[2] http://reports-archive.adm.cs.cmu.edu/anon/1995/CMU-CS-95-113.pdf

http://en.wikipedia.org/wiki/Propositional_formula
http://en.wikipedia.org/wiki/Inference_rule
http://en.wikipedia.org/wiki/List_of_rules_of_inference
"""
from collections import defaultdict

from logic import fuzzy_not, Logic, And, Or, Not

def _base_fact(atom):
    """Return the literal fact of an atom.

    Effectively, this merely strips the Not around a fact.
    """
    if isinstance(atom, Not):
        return atom.arg
    else:
        return atom

def _as_pair(atom):
    if isinstance(atom, Not):
        return (atom.arg, False)
    else:
        return (atom, True)

# XXX this prepares forward-chaining rules for alpha-network
def deduce_alpha_implications(implications):
    """deduce all implications

       Description by example
       ----------------------

       given set of logic rules:

         a -> b
         b -> c

       we deduce all possible rules:

         a -> b, c
         b -> c


       implications: [] of (a,b)
       return:       {} of a -> set([b, c, ...])
    """
    res = defaultdict(set)
    for a, b in implications:
        if a == b:
            continue    # skip a->a cyclic input

        res[a].add(b)

        # (x >> a) & (a >> b) => x >> b
        for fact in res:
            implied = res[fact]
            if a in implied:
                implied.add(b)

        # (a >> b) & (b >> x) => a >> x
        if b in res:
            res[a] |= res[b]

    # Clean up tautologies and check consistency
    for a, impl in res.iteritems():
        impl.discard(a)
        na = Not(a)
        if na in impl:
            raise ValueError('implications are inconsistent: %s -> %s %s' % (a, na, impl))

    return res


def apply_beta_to_alpha_route(alpha_implications, beta_rules):
    """apply additional beta-rules (And conditions) to already-built alpha implication tables

       TODO: write about

       - static extension of alpha-chains
       - attaching refs to beta-nodes to alpha chains


       e.g.

       alpha_implications:

       a  ->  [b, !c, d]
       b  ->  [d]
       ...


       beta_rules:

       &(b,d) -> e


       then we'll extend a's rule to the following

       a  ->  [b, !c, d, e]
    """
    x_impl = {}
    for x in alpha_implications.keys():
        x_impl[x] = (set(alpha_implications[x]), [])
    for bcond, bimpl in beta_rules:
        for bk in bcond.args:
            if bk in x_impl:
                continue
            x_impl[bk] = (set(), [])

    # static extensions to alpha rules:
    # A: x -> a,b   B: &(a,b) -> c  ==>  A: x -> a,b,c
    seen_static_extension = True
    while seen_static_extension:
        seen_static_extension = False

        for bcond, bimpl in beta_rules:
            assert isinstance(bcond, And)
            bargs = set(bcond.args)
            for x, (ximpls, bb) in x_impl.iteritems():
                x_all = ximpls | set([x])
                # A: ... -> a   B: &(...) -> a  is non-informative
                if bimpl not in x_all and bargs.issubset(x_all):
                    ximpls.add(bimpl)

                    # we introduced new implication - now we have to restore
                    # completness of the whole set.
                    bimpl_impl = x_impl.get(bimpl)
                    if bimpl_impl is not None:
                        ximpls |= bimpl_impl[0]
                    seen_static_extension=True

    # attach beta-nodes which can be possibly triggered by an alpha-chain
    for bidx, (bcond,bimpl) in enumerate(beta_rules):
        bargs = set(bcond.args)
        for x, (ximpls, bb) in x_impl.iteritems():
            x_all = ximpls | set([x])
            # A: ... -> a   B: &(...) -> a      (non-informative)
            if bimpl in x_all:
                continue
            # A: x -> a...  B: &(!a,...) -> ... (will never trigger)
            # A: x -> a...  B: &(...) -> !a     (will never trigger)
            if any(Not(xi) in bargs or Not(xi) == bimpl for xi in x_all):
                continue

            if bargs & x_all:
                bb.append(bidx)

    return x_impl


def rules_2prereq(rules):
    """build prerequisites table from rules

       Description by example
       ----------------------

       given set of logic rules:

         a -> b, c
         b -> c

       we build prerequisites (from what points something can be deduced):

         b <- a
         c <- a, b

       rules:   {} of a -> [b, c, ...]
       return:  {} of c <- [a, b, ...]

       Note however, that this prerequisites may be *not* enough to prove a
       fact. An example is 'a -> b' rule, where prereq(a) is b, and prereq(b)
       is a. That's because a=T -> b=T, and b=F -> a=F, but a=F -> b=?
    """
    prereq = defaultdict(set)
    for (a, _), impl in rules.iteritems():
        for (i, _) in impl:
            prereq[i].add(a)
    return prereq

def split_rules(rules):
    """split alpha-rules into T->T & T->F & F->T & F->F chains

       and also rewrite them to be free of not-names

       Examples
       -------

       'a' -> ['b', '!c']

       will be split into

       'a' -> ['b'] # tt: a -> b
       'a' -> ['c'] # tf: a -> !c

       and
       '!a' -> ['b']

       will become

       'b' -> ['a'] # ft: !b -> a
    """
    rel = defaultdict(set)
    for k, impl in rules.iteritems():
        if type(k) is not Not:
            for i in impl:
                rel[(k, True)].add(_as_pair(i))
                if not isinstance(i, Not):
                    rel[(i, False)].add((k, False)) # FF is related to TT
        else:
            k = k.arg
            for i in impl:
                rel[_as_pair(Not(i))].add((k, True))
    return rel


################
# RULES PROVER #
################

class TautologyDetected(Exception):
    """(internal) Prover uses it for reporting detected tautology"""
    pass


class Prover(object):
    """ai - prover of logic rules

       given a set of initial rules, Prover tries to prove all possible rules
       which follow from given premises.

       As a result proved_rules are always either in one of two forms: alpha or
       beta:

       Alpha rules
       -----------

       This are rules of the form::

         a -> b & c & d & ...


       Beta rules
       ----------

       This are rules of the form::

         &(a,b,...) -> c & d & ...


       i.e. beta rules are join conditions that say that something follows when
       *several* facts are true at the same time.
    """

    def __init__(self):
        self.proved_rules = []
        self._rules_seen  = set()

    def split_alpha_beta(self):
        """split proved rules into alpha and beta chains"""
        rules_alpha = []    # a      -> b
        rules_beta  = []    # &(...) -> b
        for a,b in self.proved_rules:
            if isinstance(a, And):
                rules_beta.append((a,b))
            else:
                rules_alpha.append((a,b) )
        return rules_alpha, rules_beta

    @property
    def rules_alpha(self):
        return self.split_alpha_beta()[0]

    @property
    def rules_beta(self):
        return self.split_alpha_beta()[1]

    def process_rule(self, a, b):
        """process a -> b rule"""   # TODO write more?
        if (not a) or isinstance(b, bool):
            return
        if isinstance(a, bool):
            return
        if (a,b) in self._rules_seen:
            return
        else:
            self._rules_seen.add((a,b))

        # this is the core of processing
        try:
            self._process_rule(a, b)
        except TautologyDetected:
            pass

    def _process_rule(self, a, b):
        # right part first

        # a -> b & c    -->  a -> b  ;  a -> c
        # (?) FIXME this is only correct when b & c != null !
        if isinstance(b, And):
            for barg in b.args:
                self.process_rule(a, barg)

        # a -> b | c    -->  !b & !c -> !a
        #               -->   a & !b -> c & !b
        #               -->   a & !c -> b & !c
        #
        # NB: the last two rewrites add 1 term, so the rule *grows* in size.
        # NB: without catching terminating conditions this could continue infinitely
        elif isinstance(b, Or):
            # detect tautology first
            if not isinstance(a, Logic):    # Atom
                # tautology:  a -> a|c|...
                if a in b.args:
                    raise TautologyDetected(a,b, 'a -> a|c|...')
            self.process_rule(And(*[Not(barg) for barg in b.args]), Not(a))

            for bidx in range(len(b.args)):
                barg = b.args[bidx]
                brest= b.args[:bidx] + b.args[bidx+1:]
                self.process_rule(And(a, Not(barg)),
                                    And(b.__class__(*brest), Not(barg)))

        # left part

        # a & b -> c    -->  IRREDUCIBLE CASE -- WE STORE IT AS IS
        #                    (this will be the basis of beta-network)
        elif isinstance(a, And):
            if b in a.args:
                raise TautologyDetected(a,b, 'a & b -> a')
            self.proved_rules.append((a,b))
            # XXX NOTE at present we ignore  !c -> !a | !b

        elif isinstance(a, Or):
            if b in a.args:
                raise TautologyDetected(a,b, 'a | b -> a')
            for aarg in a.args:
                self.process_rule(aarg, b)

        else:
            # both `a` and `b` are atoms
            self.proved_rules.append((a,b))     # a  -> b
            self.proved_rules.append((Not(b), Not(a)))   # !b -> !a

########################################

class FactRules(object):
    """Rules that describe how to deduce facts in logic space

       When defined, these rules allow implications to quickly be determined for a
       set of facts. For this precomputed deduction tables are used. see
       `deduce_all_facts`   (forward-chaining)

       Also it is possible to gather prerequisites for a fact, which is tried
       to be proven.    (backward-chaining)


       Definition Syntax
       -----------------

       a -> b       -- a=T -> b=T  (and automatically b=F -> a=F)
       a -> !b      -- a=T -> b=F
       a == b       -- a -> b & b -> a
       a -> b & c   -- a=T -> b=T & c=T
       # TODO b | c


       Internals
       ---------

       {} k -> [] of implications:

         .rel_tt      k=T -> k2=T
         .rel_tf      k=T -> k2=F
         .rel_ff      k=F -> k2=F
         .rel_ft      k=F -> k2=T

       .rel_tbeta     k=T -> [] of possibly triggering # of beta-rules
       .rel_fbeta     k=F -> ------------------//---------------------

       .rels    -- {} k -> tt, tf, ff, ft   (list of implications for k)
       .prereq  -- {} k <- [] of k's prerequisites

       .defined_facts -- set of defined fact names
    """

    def __init__(self, rules):
        """Compile rules into internal lookup tables"""

        if isinstance(rules, basestring):
            rules = rules.splitlines()

        # --- parse and process rules ---
        P = Prover()

        for rule in rules:
            # XXX `a` is hardcoded to be always atom
            a, op, b = rule.split(None, 2)

            a = Logic.fromstring(a)
            b = Logic.fromstring(b)

            if op == '->':
                P.process_rule(a, b)
            elif op == '==':
                P.process_rule(a, b)
                P.process_rule(b, a)
            else:
                raise ValueError('unknown op %r' % op)

        # --- build deduction networks ---

        # deduce alpha implications
        impl_a = deduce_alpha_implications(P.rules_alpha)

        # now:
        # - apply beta rules to alpha chains  (static extension), and
        # - further associate beta rules to alpha chain (for inference at runtime)
        impl_ab = apply_beta_to_alpha_route(impl_a, P.rules_beta)

        # extract defined fact names
        self.defined_facts = set(_base_fact(k) for k in impl_ab.keys())

        # now split each rule into four logic chains
        # (removing betaidxs from impl_ab view) (XXX is this needed?)
        impl_ab_ = dict( (k,impl)  for k, (impl,betaidxs) in impl_ab.iteritems())
        rel_alpha = split_rules(impl_ab_)
        rel_beta = {}
        for k, (impl,betaidxs) in impl_ab.iteritems():
            rel_beta[_as_pair(k)] = betaidxs

        self.beta_rules = P.rules_beta

        # build rels (forward chains)
        K = set(rel_alpha) | set(rel_beta)

        rels = {}
        empty= ()
        for kv in K:
            rels[kv] = tuple([rule.get(kv, empty) for rule in rel_alpha, rel_beta])

        self.rels = rels

        # build prereq (backward chains)
        prereq = defaultdict(set)
        rel_prereq = rules_2prereq(rel_alpha)
        for k, pitems in rel_prereq.iteritems():
            prereq[k] |= pitems
        self.prereq = prereq

class FactKB(dict):
    def __init__(self, rules):
        self.rules = rules

    def copy(self):
        new = self.__class__(self.rules)
        new.update(self)
        return new

    def deduce_all_facts(self, facts):
        """Deduce all facts from known facts ({} or [] of (k,v))

           *********************************************
           * This is the workhorse, so keep it *fast*. *
           *********************************************

           base  --  previously known facts (must be: fully deduced set)
                     attention: base is modified *in place*  /optional/

           providing `base` could be needed for performance reasons -- we don't
           want to spend most of the time just re-deducing base from base
           (e.g. #base=50, #facts=2)
        """
        # keep frequently used attributes locally, so we'll avoid extra
        # attribute access overhead
        rels = self.rules.rels
        beta_rules = self.rules.beta_rules

        def x_new_facts(facts):
            for k, v in facts:
                if k in self and self[k] is not None:
                    assert self[k] == v, \
                            ('inconsistency between facts', self, k, v)
                    continue
                else:
                    self[k] = v

        if isinstance(facts, dict):
            fseq = facts.iteritems()
        else:
            fseq = facts

        while True:
            beta_maytrigger = set()

            # --- alpha chains ---
            for k, v in fseq:
                #new_fact(k, v)
                if k in self and self[k] is not None:
                    assert self[k] == v, \
                            ('inconsistency between facts', self, k, v)
                    # performance-wise it is important not to fire implied rules
                    # for already-seen fact -- we already did them all.
                    continue
                else:
                    self[k] = v

                # some known fact -- let's follow its implications
                if v is not None:
                    # lookup routing tables
                    try:
                        alpha, beta = rels[(k, v)]
                    except KeyError:
                        pass
                    else:
                        # Now we have routing tables with *all* the needed
                        # implications for this k. This means we do not have to
                        # process each implications recursively!
                        # XXX this ^^^ is true only for alpha chains
                        x_new_facts(alpha)

                        beta_maytrigger.update(beta)

            # --- beta chains ---

            # if no beta-rules may trigger -- it's an end-of-story
            if not beta_maytrigger:
                break
            fseq = []
            # let's see which beta-rules to trigger
            for bidx in beta_maytrigger:
                bcond,bimpl = beta_rules[bidx]
                # let's see whether bcond is satisfied
                for bk in bcond.args:
                    try:
                        if type(bk) is Not:
                            bv = fuzzy_not(self[bk.arg])
                        else:
                            bv = self[bk]
                    except KeyError:
                        break   # fact not found -- bcond not satisfied
                    # one of bcond's condition does not hold
                    if not bv:
                        break
                else:
                    # all of bcond's condition hold -- let's fire this beta rule
                    if type(bimpl) is Not:
                        bimpl = bimpl.arg
                        v = False
                    else:
                        v = True
                    fseq.append( (bimpl,v) )
        return self
