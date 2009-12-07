# -*- coding: utf-8 -*-

"""This is rule-based deduction system for SymPy

The whole thing is split into two parts

 - rules compilation and preparation of tables
 - runtime inference

For rule-based inference engines, the classical work is RETE algorithm [1], [2]
Although we are not implementing it in full (or even significantly) it's still
still worth a read to understand the underlying ideas.

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



from logic import fuzzy_not, name_not, Logic, And, Not

def list_populate(l, item, skipif=None):
    """update list with an item, but only if it is not already there"""
    if item != skipif and (item not in l):
        l.append(item)



# XXX this prepares forward-chaining rules for alpha-network

# XXX better words?
# NB: this procedure does not 'know' about not-relations
# i.e. it treats 'a' and '!a' as independent variables while deriving implications.
#
# only in the end, we check that the whole result is consitent
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
       return:       {} of a -> [b, c, ...]
    """

    res = {} # a -> [] of implications(a)

    # NOTE: this could be optimized with e.g. toposort (?), but we need this
    # only at FactRules creation time (i.e. only once), so there is no demand
    # in further optimising this.
    for a,b in implications:
        if a==b:
            continue    # skip a->a cyclic input

        I = res.setdefault(a,[])
        list_populate(I,b)

        # UC:  -------------------------
        #     |                         |
        #     v                         |
        # a='rat' -> b='real' ==> (a_='int') -> 'real'
        for a_ in res:
            ra_ = res[a_]
            if a in ra_:
                list_populate(ra_, b, skipif=a_)


        # UC:
        # a='pos' -> b='real' && (already have b='real' -> 'complex')
        #                   ||
        #                   vv
        # a='pos' -> 'complex'
        if b in res:
            ra = res[a]
            for b_ in res[b]:
                list_populate(ra, b_, skipif=a)

    #print 'D:', res

    # let's see if the result is consistent
    for a, impl in res.iteritems():
        na = name_not(a)
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
        x_impl[x] = (alpha_implications[x][:], [])

    # let's ensure that every beta-rule has appropriate entry (maybe even
    # empty) in the table.
    for bidx, (bcond,bimpl) in enumerate(beta_rules):
        for bk in bcond.args:
            if bk in x_impl:
                continue

            x_impl[bk] = ([], [])


    # XXX maybe simplify passes logic?

    # we do it in 2 phases:
    #
    # 1st phase -- only do static extensions to alpha rules
    # 2nd phase -- attach beta-nodes which can be possibly triggered by an
    #              alpha-chain
    phase=1
    while True:
        seen_static_extension=False

        for bidx, (bcond,bimpl) in enumerate(beta_rules):
            assert isinstance(bcond, And)

            for x, (ximpls, bb) in x_impl.iteritems():
                # A: x -> a     B: &(...) -> x      (non-informative)
                if x == bimpl:  # XXX bimpl may become a list
                    continue

                # A: ... -> a   B: &(...) -> a      (non-informative)
                if bimpl in ximpls:
                    continue

                # A: x -> a     B: &(a,!x) -> ...   (will never trigger)
                if Not(x) in bcond.args:
                    continue

                # A: x -> a...  B: &(!a,...) -> ... (will never trigger)
                # A: x -> a...  B: &(...) -> !a     (will never trigger)
                maytrigger=True
                for xi in ximpls:
                    if (Not(xi) in bcond.args)  or  (Not(xi) == bimpl):
                        maytrigger=False
                        break

                if not maytrigger:
                    continue

                # A: x -> a,b   B: &(a,b) -> c      (static extension)
                #                                       |
                # A: x -> a,b,c <-----------------------+
                for barg in bcond.args:
                    if not ( (barg == x) or (barg in ximpls) ):
                        break
                else:
                    assert phase==1
                    list_populate(ximpls, bimpl)    # XXX bimpl may become a list

                    # we introduced new implication - now we have to restore
                    # completness of the whole set.
                    bimpl_impl = x_impl.get(bimpl)
                    if bimpl_impl is not None:
                        for _ in bimpl_impl[0]:
                            list_populate(ximpls, _)

                    seen_static_extension=True
                    continue


                # does this beta-rule even has a chance to be triggered ?
                if phase == 2:
                    for barg in bcond.args:
                        if (barg == x) or (barg in ximpls):
                            bb.append( bidx )
                            break

        # no static extensions was seen at this pass -- lets move to phase2
        if phase==1 and (not seen_static_extension):
            phase = 2
            continue

        # let's finish at the end of phase2
        if phase==2:
            break


    return x_impl



# XXX this is something related to backward-chaining for alpha-network (?)
# TODO adapt to handle beta nodes ?

# TODO we need to sort prerequisites somehow, so that more closer prerequisites
# would be tried first. e.g.
#
# nonzero <- [zero, integer, ..., prime], but not
# nonzero <- [prime, ...]
#
# this could be done in other place...
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

    prereq = {}
    for a, impl in rules.iteritems():
        for i in impl:
            pa = prereq.setdefault(i,[])
            pa.append(a)

    return prereq



def split_rules_tt_tf_ft_ff(rules):
    """split alpha-rules into T->T & T->F & F->T & F->F chains

       and also rewrite them to be free of not-names

       Example
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

#   print 'split_rules_tt_tf_ft_ff'
#   print rules

    tt = {}
    tf = {}
    ft = {}

    for k,impl in rules.iteritems():
        # k is not not-name
        if k[:1] != '!':
            for i in impl:
                if i[:1] != '!':
                    dd = tt
                else:
                    dd = tf
                    i  = i[1:]

                I = dd.setdefault(k,[])
                list_populate(I, i)

        # k is not-name
        else:
            k = k[1:]

            for i in impl:
                if i[:1] != '!':
                    dd = ft
                else:
                    dd = tt
                    i  = i[1:]

                I = dd.setdefault(i,[])
                list_populate(I, k)


    # FF is related to TT
    ff = {}
    for k,impl in tt.iteritems():
        for i in impl:
            I = ff.setdefault(i,[])
            I.append(k)

    return tt, tf, ft, ff





################
# RULES PROVER #
################

# XXX only for debugging -- kill me ?
dbg_level = 0


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

    __slots__ = ['proved_rules',    # [] of a,b     (a -> b  rule)
                 '_rules_seen',     # setof all seen rules
                ]


    def __init__(self):
        self.proved_rules = []
        self._rules_seen  = set()

    def print_proved(self, title='proved rules'):
        print '\n--- %s ---' % title
        for a,b in self.proved_rules:
            print '%s\t->  %s' % (a,b)

        print '   - - - - -   '
        print


    def print_beta(self, title='proved rules (beta)'):
        print '\n --- %s ---' % title

        for n, (a,b) in enumerate(self.rules_beta):
            print '[#%i]  %s\t->  %s' % (n,a,b)

        print '   - - - - -   '
        print


    def split_alpha_beta(self):
        """split proved rules into alpha and beta chains"""
        rules_alpha = []    # a      -> b
        rules_beta  = []    # &(...) -> b

        for a,b in self.proved_rules:
            if isinstance(a, And):
                rules_beta.append( (a,b) )
            else:
                rules_alpha.append((a,b) )

        return rules_alpha, rules_beta


    @property
    def rules_alpha(self):
        return self.split_alpha_beta()[0]

    @property
    def rules_beta(self):
        return self.split_alpha_beta()[1]



    # XXX rename?
    def process_rule(self, a, b):
        """process a -> b rule"""   # TODO write more?

        # XXX this seems to be needed, but I can't get it to work without cycles
        # XXX sigh ...
        #if isinstance(a, And):
        #    a = a.expand()
        #if isinstance(b, And):
        #    b = b.expand()

#       print ' -- %s\t->  %s' % (a,b)

        # we don't want rules whose premise is F, or who are trivial
        # skip:  F -> ...
        # skip:  ... -> F
        # skip:  ... -> T
        if (not a) or isinstance(b, bool):
            return

        # XXX can't we perform without this?
        # skip:  T -> ...
        if isinstance(a, bool):
            return


        # cycle detection: we don't want to process already seen rules
        if (a,b) in self._rules_seen:
            # already seen rule
            # XXX is it possible to design process_rule so we don't need
            # XXX this recursion detector? For me current answer is 'no'
            #                                                   -- kirr
            return
        else:
            self._rules_seen.add((a,b))

        # this is the core of processing
        try:
            self._process_rule(a, b)
        except TautologyDetected, t:
            #print 'Tautology:  %s -> %s   (%s)' % (t.args[0], t.args[1], t.args[2])
            pass

        else:
            # XXX negating here seems like the way to go. However we do
            # XXX *not* do it -- read below about why ...
            #
            # if not tautology, process negative rule as well
            # a -> b  -->  !b -> !a
        #   self.process_rule(Not(b), Not(a))


            # Unfortunately we can't handle the above negation -- in the
            # end our cycles detector is unable to catch all cycles,
            # especially ones where noise is mixed to already-known-facts,
            # e.g.
            #
            # say
            #
            #   a -> b
            #
            # is known fact, and at some point we get
            #
            #   a -> |(....., b)
            #
            # At present we are unable to detect this situation.
            #
            # I'm sorry -- I have no deep understanding of this logic
            # processing stuff and I'm very short on time ...
            #                                                  -- kirr
            pass


    # --- processing ---

    def _process_rule(self, a, b):

        # right part first
        if isinstance(b, Logic):
            # a -> b & c    -->  a -> b  ;  a -> c
            # (?) FIXME this is only correct when b & c != null !
            if   b.op == '&':
                for barg in b.args:
                    self.process_rule(a, barg)

            # a -> b | c    -->  !b & !c -> !a
            #               -->   a & !b -> c & !b
            #               -->   a & !c -> b & !c
            #
            # NB: the last two rewrites add 1 term, so the rule *grows* in size.
            # NB: without catching terminating conditions this could continue infinitely
            elif b.op == '|':
                # detect tautology first
                if not isinstance(a, Logic):    # Atom
                    # tautology:  a -> a|c|...
                    if a in b.args:
                        raise TautologyDetected(a,b, 'a -> a|c|...')

                self.process_rule(
                    And( *[Not(barg) for barg in b.args] ),
                    Not(a)
                )

                for bidx in range(len(b.args)):
                    barg = b.args[bidx]
                    brest= b.args[:bidx] + b.args[bidx+1:]

                    self.process_rule(
                        And( a, Not(barg) ),
                        And( b.__class__(*brest), Not(barg) )
                    )

            else:
                raise ValueError('unknown b.op %r' % b.op)


        # left part
        elif isinstance(a, Logic):
            # a & b -> c    -->  IRREDUCIBLE CASE -- WE STORE IT AS IS
            #                    (this will be the basis of beta-network)
            if   a.op == '&':
                # at this stage we should have right part alredy in simple form
                assert not isinstance(b, Logic)

                # tautology:  a & b -> a
                if b in a.args:
                    raise TautologyDetected(a,b, 'a & b -> a')

                self.proved_rules.append((a,b))
                # XXX NOTE at present we ignore  !c -> !a | !b

            # a | b -> c    -->  a -> c  ;  b -> c
            elif a.op == '|':

                # tautology:  a | b -> a
                if b in a.args:
                    raise TautologyDetected(a,b, 'a | b -> a')

                for aarg in a.args:
                    self.process_rule(aarg, b)

            else:
                raise ValueError('unknown a.op %r' % a.op)

        else:
            # now this is the case where both `a` and `b` are atoms
            na, nb = name_not(a), name_not(b)

            self.proved_rules.append((a,b))     # a  -> b

            # XXX one day we may not need this -- read comments about
            # XXX negation close to the end of `process_rule`
            self.proved_rules.append((nb,na))   # !b -> !a


    def dbg_process_rule_2(a, b):
        global dbg_level
        print '%s%s\t->  %s' % (' '*(2*dbg_level), a, b)
        dbg_level += 1

        try:
            old_process_rule_2(a, b)
        finally:
            dbg_level -= 1

    #old_process_rule_2 = process_rule_2
    #process_rule_2 = dbg_process_rule_2



########################################


# TODO link to RETE ?
# TODO link to proposition logic ?
class FactRules:
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

            if   op == '->':
                P.process_rule(a, b)
            elif op == '==':
                P.process_rule(a, b)
                P.process_rule(b, a)
            else:
                raise ValueError('unknown op %r' % op)


        #P.print_proved('RULES')
        #P.print_beta('BETA-RULES')


        # --- build deduction networks ---

        # deduce alpha implications
        impl_a = deduce_alpha_implications(P.rules_alpha)

        # now:
        # - apply beta rules to alpha chains  (static extension), and
        # - further associate beta rules to alpha chain (for inference at runtime)
        impl_ab = apply_beta_to_alpha_route(impl_a, P.rules_beta)

        if 0:
            print '\n --- ALPHA-CHAINS (I) ---'

            for a,b in impl_a.iteritems():
                print '%s\t->  α(%2i):%s' % (a,len(b),b)

            print '   - - - - -   '
            print

        if 0:
            print '\n --- ALPHA-CHAINS (II) ---'

            for a,(b,bb) in impl_ab.iteritems():
                print '%s\t->  α(%2i):%s  β(%s)' % (a,len(b),b, ' '.join(str(x) for x in bb))

            print '   - - - - -   '
            print


        # extract defined fact names
        self.defined_facts = set()

        for k in impl_ab.keys():
            if k[:1] == '!':
                k = k[1:]

            self.defined_facts.add(k)

        #print 'defined facts:  (%2i)  %s' % (len(self.defined_facts), self.defined_facts)


        # now split each rule into four logic chains
        # (removing betaidxs from impl_ab view) (XXX is this needed?)
        impl_ab_ = dict( (k,impl)  for k, (impl,betaidxs) in impl_ab.iteritems())
        rel_tt, rel_tf, rel_ft, rel_ff = split_rules_tt_tf_ft_ff(impl_ab_)

        # XXX merge me with split_rules_tt_tf_ft_ff ?
        rel_tbeta = {}
        rel_fbeta = {}

        for k, (impl,betaidxs) in impl_ab.iteritems():
            if k[:1] == '!':
                rel_xbeta = rel_fbeta
                k         = name_not(k)
            else:
                rel_xbeta = rel_tbeta

            rel_xbeta[k] = betaidxs



        self.rel_tt = rel_tt
        self.rel_tf = rel_tf
        self.rel_tbeta  = rel_tbeta
        self.rel_ff = rel_ff
        self.rel_ft = rel_ft
        self.rel_fbeta  = rel_fbeta

        self.beta_rules = P.rules_beta

        # build rels (forward chains)
        K = set (rel_tt.keys())
        K.update(rel_tf.keys())
        K.update(rel_ff.keys())
        K.update(rel_ft.keys())

        rels = {}
        empty= ()
        for k in K:
            tt = rel_tt.get(k,empty)
            tf = rel_tf.get(k,empty)
            ft = rel_ft.get(k,empty)
            ff = rel_ff.get(k,empty)

            tbeta = rel_tbeta.get(k,empty)
            fbeta = rel_fbeta.get(k,empty)

            rels[k] = tt, tf, tbeta,  ft, ff, fbeta


        self.rels = rels

        # build prereq (backward chains)
        prereq = {}
        for rel in [rel_tt, rel_tf, rel_ff, rel_ft]:
            rel_prereq = rules_2prereq(rel)
            for k,pitems in rel_prereq.iteritems():
                kp = prereq.setdefault(k,[])
                for p in pitems:
                    list_populate(kp, p)

        self.prereq = prereq



    # --- DEDUCTION ENGINE: RUNTIME CORE ---

    # TODO: add proper support for U (None), i.e.
    #       integer=U  -> rational=U  ??? (XXX i'm not sure)


    def deduce_all_facts(self, facts, base=None):
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
        rels = self.rels
        beta_rules = self.beta_rules

        if base is not None:
            new_facts = base
        else:
            new_facts = {}

        # XXX better name ?
        def x_new_facts(keys, v):
            for k in keys:
                if k in new_facts and new_facts[k] is not None:
                    assert new_facts[k] == v, \
                            ('inconsistency between facts', new_facts, k, v)
                    continue
                else:
                    new_facts[k] = v

        if type(facts) is dict:
            fseq = facts.iteritems()
        else:
            fseq = facts


        while True:
            beta_maytrigger = set()

            # --- alpha chains ---
            #print '**'
            for k,v in fseq:
                #print '--'

                # first, convert name to be not a not-name
                if k[:1] == '!':
                    k = name_not(k)
                    v = fuzzy_not(v)

                #new_fact(k, v)
                if k in new_facts:
                    assert new_facts[k] == v, \
                            ('inconsistency between facts', new_facts, k, v)
                    # performance-wise it is important not to fire implied rules
                    # for already-seen fact -- we already did them all.
                    continue
                else:
                    new_facts[k] = v

                # some known fact -- let's follow its implications
                if v is not None:
                    # lookup routing tables
                    try:
                        tt, tf, tbeta,  ft, ff, fbeta = rels[k]
                    except KeyError:
                        pass
                    else:
                        # Now we have routing tables with *all* the needed
                        # implications for this k. This means we do not have to
                        # process each implications recursively!
                        # XXX this ^^^ is true only for alpha chains

                        # k=T
                        if v:
                            x_new_facts(tt, True)   # k -> i
                            x_new_facts(tf, False)  # k -> !i

                            beta_maytrigger.update(tbeta)

                        # k=F
                        else:
                            x_new_facts(ft, True)   # !k -> i
                            x_new_facts(ff, False)  # !k -> !i

                            beta_maytrigger.update(fbeta)


            # --- beta chains ---

            # if no beta-rules may trigger -- it's an end-of-story
            if not beta_maytrigger:
                break
            #print '(β) MayTrigger: %s' % beta_maytrigger

            fseq = []

            # XXX this is dumb proof-of-concept trigger -- we'll need to optimize it
            # let's see which beta-rules to trigger
            for bidx in beta_maytrigger:
                bcond,bimpl = beta_rules[bidx]

                # let's see whether bcond is satisfied
                for bk in bcond.args:
                    try:
                        if bk[:1] == '!':
                            bv = fuzzy_not(new_facts[bk[1:]])
                        else:
                            bv = new_facts[bk]
                    except KeyError:
                        break   # fact not found -- bcond not satisfied

                    # one of bcond's condition does not hold
                    if not bv:
                        break

                else:
                    # all of bcond's condition hold -- let's fire this beta rule
                    #print '(β) Trigger #%i  (%s)' % (bidx, bimpl)

                    if bimpl[:1] == '!':
                        bimpl = bimpl[1:]
                        v = False
                    else:
                        v = True
                    fseq.append( (bimpl,v) )
        return new_facts


########################################



# this was neeeded for testing -- feel free to remove this when everything
# settles.
#                                                               -- kirr
if 0:

    FF = FactRules([

        'integer        ->  rational',
        'rational       ->  real',
        'real           ->  complex',
        'imaginary      ->  complex',
        'complex        ->  commutative',

        'odd            ==  integer & !even',
        'even           ==  integer & !odd',

        'real           ==  negative | zero | positive',

        'positive       ->  real & !negative & !zero',
        'negative       ->  real & !positive & !zero',

        'nonpositive    ==  real & !positive',
        'nonnegative    ==  real & !negative',

        'zero           ->  infinitesimal & even',

        'prime          ->  integer & positive',
        'composite      ==  integer & positive & !prime',

        'irrational     ==  real & !rational',

        'imaginary      ->  !real',     # XXX new ?


        '!bounded     ==  unbounded',
        '!commutative ==  noncommutative',
        '!complex     ==  noncomplex',
    #   '!integer     ==  noninteger',
        'noninteger     ==  real & !integer',   # XXX ok?
        '!zero        ==  nonzero',
        '!homogeneous ==  inhomogeneous',

        # XXX do we need this ?
        'finite     ->  bounded',       # XXX do we need this?
        'finite     ->  !zero',         # XXX wrong?
        'infinitesimal ->  !finite',    # XXX is this ok?

        # TODO we should remove this (very seldomly used, but affect performance):
        'nni    ==  integer & nonnegative',
        'npi    ==  integer & nonpositive',
        'pi     ==  integer & positive',
        'ni     ==  integer & negative',
    ])

    """
    FF = FactRules([
        'real   == pos | npos',

        'pos    -> real & !npos',
        'npos   -> real & !pos'
    ])
    """

    """
    ff = FactRules([
        'real   ->  complex',
        'real   ==  neg | zero | pos',

        'neg    ->  real & !pos & !zero',
        'zero   ->  real & !neg & !pos',
        'pos    ->  real & !neg & !zero',

        'npos   ==  real & !pos',
        'real   ==  pos | npos' # XXX this is a workaround for bug in process_rule
    ])

    FF = ff
    """
