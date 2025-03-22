### **fix**(lambdify): make dummify=None default, replace only invalid identifiers

<!-- Your title above should be a short description of what
was changed. Do not include the issue number in the title. -->

### **References** to other Issues or PRs

<!-- If this pull request fixes an issue, write "Fixes #NNNN" in that exact
format, e.g. "Fixes #1234" (see
https://tinyurl.com/auto-closing for more information). Also, please
write a comment on that issue linking back to this pull request once it is
open. -->

Fixes #12463

### **Description**  
Fixes #12463  
Currently, `lambdify` uses `dummify=True` by default, which replaces all symbols with `_Dummy_xxx`, breaking function introspection. This PR changes the default to `dummify=None`, ensuring that only invalid identifiers or Python keywords are replaced.

### **Technical Changes**  
1. Changed the default value of `dummify` from `True` to `None` in `lambdify`.  
2. Added identifier validation logic: symbols are replaced only if their names are invalid Python identifiers or reserved keywords.  
3. Updated test cases to cover valid and invalid identifier scenarios.  

### **Testing**  
- Passed local tests with `./bin/test sympy/utilities`.  
- Added `test_dummify_smart` to verify the smart replacement logic.  

### **Compatibility**  
- Backward compatible: Users can still force replacement with `dummify=True` or disable it with `dummify=False`.  

### **Related Issues**  
Closes #12463  

<!-- BEGIN RELEASE NOTES -->

### **Description**  
Fixes #12463  
Currently, `lambdify` uses `dummify=True` by default, which replaces all symbols with `_Dummy_xxx`, breaking function introspection. This PR changes the default to `dummify=None`, ensuring that only invalid identifiers or Python keywords are replaced.

### **Technical Changes**  
1. Changed the default value of `dummify` from `True` to `None` in `lambdify`.  
2. Added identifier validation logic: symbols are replaced only if their names are invalid Python identifiers or reserved keywords.  
3. Updated test cases to cover valid and invalid identifier scenarios.  

### **Testing**  
- Passed local tests with `./bin/test sympy/utilities`.  
- Added `test_dummify_smart` to verify the smart replacement logic.  

### **Compatibility**  
- Backward compatible: Users can still force replacement with `dummify=True` or disable it with `dummify=False`.  

### **Related Issues**  
Closes #12463  

<!-- BEGIN RELEASE NOTES -->

### **Description**  
Fixes #12463  
Currently, `lambdify` uses `dummify=True` by default, which replaces all symbols with `_Dummy_xxx`, breaking function introspection. This PR changes the default to `dummify=None`, ensuring that only invalid identifiers or Python keywords are replaced.

### **Technical Changes**  
1. Changed the default value of `dummify` from `True` to `None` in `lambdify`.  
2. Added identifier validation logic: symbols are replaced only if their names are invalid Python identifiers or reserved keywords.  
3. Updated test cases to cover valid and invalid identifier scenarios.  

### **Testing**  
- Passed local tests with `./bin/test sympy/utilities`.  
- Added `test_dummify_smart` to verify the smart replacement logic.  

### **Compatibility**  
- Backward compatible: Users can still force replacement with `dummify=True` or disable it with `dummify=False`.  

### **Related Issues**  
Closes #12463  

<!-- BEGIN RELEASE NOTES -->

- Fixed a bug in `lambdify` where `dummify=True` broke function introspection.

- Changed the default value of `dummify` to `None` to replace only invalid identifiers.

- Added tests for the new `dummify=None` behavior.

<!-- END RELEASE NOTES -->
