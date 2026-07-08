# Instructions for AI agents

## Opening pull requests

When opening a PR for this repo (sympy/sympy), always account for AI usage as
required by `.github/PULL_REQUEST_TEMPLATE.md` and sympy's
[Policy on AI Generated Code](doc/src/contributing/ai-generated-code-policy.rst):

- **Always fill in the `#### AI Generation Disclosure` section** of the PR
  template. Disclose the tool or agent used (and the model) and specify
  which files or line ranges were AI-generated. Disclosure is not required for
  minor assistive tasks such as spell-checking or code review in primarily
  human-authored work — but if in doubt, disclose. Only write "NO AI USE" when
  there was genuinely no AI-generated code or text.
- **Never leave that section blank and never delete it**, or the PR will be
  closed automatically.
- **Do not use AI-generated text for the PR description prose** (the "Brief
  description", "Other comments", etc.). The policy forbids it and the PR will
  be closed. Write those sections in your own words as a plain factual summary,
  not model-generated prose.
- **Keep every section heading** from the template intact (References, Brief
  description, Other comments, AI Generation Disclosure, Release Notes).
