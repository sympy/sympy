from sympy import Matrix
#from sympy.utilities.pytest import raises

from sympy.tensor import *
from sympy import sin, symbols, cos

def test_df_varlist():
    x1, x2, x3 = symbols('x1 x2 x3')
    f=x1**2*x2+sin(x2*x3-x2)
        
    #var-list
    var_list=[x1, x2, x3]
    
    print ('test_df_varlist  <============================ actual test code')
    assert df(f, var_list) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    print ('test_df_varlist_type  <============================ actual test code')
    assert isinstance(df(f, var_list), list)    

    print ('test_df_varlist_l  <============================ actual test code')
    assert df(f, var_list, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    print ('test_df_varlist_l_type  <============================ actual test code')
    assert isinstance(df(f, var_list, 'l'), list)

    print ('test_df_varlist_a  <============================ actual test code')
    assert df(f, var_list, 'a') == list2arraypy([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_varlist_a_type  <============================ actual test code')
    assert isinstance(df(f, var_list, 'a'), Arraypy)
    
    print ('test_df_varlist_t  <============================ actual test code')
    assert df(f, var_list, 't') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]) 
    assert isinstance(df(f, var_list, 't'),Tensor)
    assert df(f, var_list, 't').type_pq == (0,1)
    
def test_df_var_tnsr0():
    x1, x2, x3 = symbols('x1 x2 x3')
    f=x1**2*x2+sin(x2*x3-x2)            
                
    #var-tensor 0
    var_tnsr0=Tensor(Arraypy(3),(1))
    var_tnsr0[0]=x1
    var_tnsr0[1]=x2
    var_tnsr0[2]=x3 
    
    print ('test_df_var_tnsr0  <============================ actual test code')
    assert df(f, var_tnsr0) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(df(f, var_tnsr0), list)     
        
    print ('test_df_var_tnsr0_l  <============================ actual test code')
    assert df(f, var_tnsr0, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(df(f, var_tnsr0, 'l'), list)    
    
    print ('test_df_var_tnsr0_a  <============================ actual test code')
    assert df(f, var_tnsr0, 'a') == list2arraypy([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_var_tnsr0_a_type  <============================ actual test code')
    assert isinstance(df(f, var_tnsr0, 'a'), Arraypy)
    
    print ('test_df_var_tnsr0_t  <============================ actual test code')
    assert df(f, var_tnsr0, 't') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    print ('test_df_var_tnsr0_t_type  <============================ actual test code')
    assert isinstance(df(f, var_tnsr0, 't'), Tensor)
    assert df(f, var_tnsr0, 't').type_pq == (0,1)  
      
def test_df_var_tnsr1():
        x1, x2, x3 = symbols('x1 x2 x3')
        f=x1**2*x2+sin(x2*x3-x2)            
                    
        #var-tensor 1
        var_tnsr1=Arraypy([1,3,1]).to_tensor(1)
        var_tnsr1[1]=x1
        var_tnsr1[2]=x2
        var_tnsr1[3]=x3 
        
        print ('test_df_var_tnsr1  <============================ actual test code')
        assert df(f, var_tnsr1) == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
        assert isinstance(df(f, var_tnsr1), list)        
            
        print ('test_df_var_tnsr1_l  <============================ actual test code')
        assert df(f, var_tnsr1, 'l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
        assert isinstance(df(f, var_tnsr1, 'l'), list)    
        
        """
        print ('test_df_var_tnsr1_a  <============================ actual test code')
        assert df(f, var_tnsr1, 'a') == list2arraypy([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
        
        print ('test_df_var_tnsr1_a_type  <============================ actual test code')
        assert isinstance(df(f, var_tnsr1, 'a'), Arraypy)
        
        print ('test_df_var_tnsr1_t  <============================ actual test code')
        assert df(f, var_tnsr1, 't') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
        print ('test_df_var_tnsr1_t_type  <============================ actual test code')
        assert isinstance(df(f, var_tnsr1, 't'), Tensor)
        assert df(f, var_tnsr1, 't').type_pq == (0,1)
        """
        
def test_div_var_x_list():
    
    x1, x2, x3 = symbols('x1 x2 x3')
    
    X = [x1*x2**3,x2-cos(x3),x3**3-x1]
    var= [x1, x2, x3]
    g=Matrix([[2,1,0],[1,3,0],[0,0,1]])
    
    ten=Arraypy([2,3,0]).to_tensor((-1,-1))
    ten[0,0]=2
    ten[0,1]=1
    ten[0,2]=0
    ten[1,0]=1
    ten[1,1]=3
    ten[1,2]=0
    ten[2,0]=0
    ten[2,1]=0
    ten[2,2]=1
    
    #g задано tensor, индекс с 1 и var-list
    ten1=Arraypy([2,3,1]).to_tensor((-1,-1))
    ten1[1,1]=2
    ten1[1,2]=1
    ten1[1,3]=0
    ten1[2,1]=1
    ten1[2,2]=3
    ten1[2,3]=0
    ten1[3,1]=0
    ten1[3,2]=0
    ten1[3,3]=1 
    
    print ('test_div_var_x_list  <============================ actual test code')
    assert div(X, var) == x2**3 + 3*x3**2 + 1     
    
    print ('test_div_var_x_list_g  <============================ actual test code')
    assert div(X, var, g) == x2**3 + 3*x3**2 + 1     
    
    print ('test_div_var_x_list_t0  <============================ actual test code')
    assert div(X, var, ten) == x2**3 + 3*x3**2 + 1  
    
    print ('test_div_var_x_list_t1  <============================ actual test code')
    assert div(X, var, ten1) == x2**3 + 3*x3**2 + 1 
    
def test_grad():
    x1, x2, x3= symbols('x1 x2 x3')
    f=x1**2*x2+sin(x2*x3-x2)
    var1=[x1, x2, x3] 
    
    k1=Arraypy([1,3,0]).to_tensor(1)
    k1[0]=x1
    k1[1]=x2
    k1[2]=x3   
    
    #g задано tensor, индекс с 1 и var-list
    a=Arraypy([2,3,1])
    b=a.to_tensor((-1,-1))
    b[1,1]=2
    b[1,2]=1
    b[1,3]=0
    b[2,1]=1
    b[2,2]=3
    b[2,3]=0
    b[3,1]=0
    b[3,2]=0
    b[3,3]=1 
    
    print ('test_grad  <============================ actual test code')
    assert grad(f,var1,output_type='l') == [2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)]
    assert isinstance(grad(f,var1,output_type='l'), list) 
    
    """
    print ('test_grad  <============================ actual test code')
    assert grad(f,var1,output_type='t') == list2tensor([2*x1*x2, x1**2 + (x3 - 1)*cos(x2*x3 - x2), x2*cos(x2*x3 - x2)])
    assert isinstance(grad(f,var1,output_type='t'), Tensor)    
    assert grad(f,var1,output_type='t').type_pq == (1,0)
    """
    
    
   
    
   
            
