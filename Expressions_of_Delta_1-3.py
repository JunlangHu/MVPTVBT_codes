# Ultimate interactive calculator
import sympy as sp

def get_D_expressions(m1_expr_str, s_val_str=None):
    """Get D expressions based on the input m1 expression"""
    try:
        # Define symbols
        m1_sym, sigma_sym, s_sym = sp.symbols('m_1 sigma s')
        
        # Parse m1 expression
        m1_expr = sp.sympify(m1_expr_str, locals={'sigma': sigma_sym, 's': s_sym})
        
        # Parse s value (if provided)
        s_val = sp.sympify(s_val_str) if s_val_str else s_sym
        
        # D1 expression
        D1 = sigma_sym**2 / 2
        
        # D2 expression (descending powers of m1)
        D2_coeffs = [
            -1,
            2*s_sym,
            -s_sym**2 - 3*sigma_sym**2/2 + 2,
            2*s_sym*sigma_sym**2 - 2*s_sym,
            -s_sym**2*sigma_sym**2/2 - 3*sigma_sym**4/4 + 3*sigma_sym**2/2 - 1,
            s_sym*sigma_sym**4/2 - s_sym*sigma_sym**2/2,
            -sigma_sym**6/8 + sigma_sym**4/4
        ]
        D2 = sum(coeff * m1_sym**(6-i) for i, coeff in enumerate(D2_coeffs))
        
        # D3 expression (descending powers of m1)
        D3_coeffs = [
            -3*sigma_sym**2/2,
            4*s_sym*sigma_sym**2,
            -7*s_sym**2*sigma_sym**2/2 - 3*sigma_sym**4 + 3*sigma_sym**2,
            s_sym**3*sigma_sym**2 + 6*s_sym*sigma_sym**4 - 5*s_sym*sigma_sym**2,
            -7*s_sym**2*sigma_sym**4/2 + 2*s_sym**2*sigma_sym**2 - 9*sigma_sym**6/4 + 3*sigma_sym**4 - 3*sigma_sym**2/2,
            s_sym**3*sigma_sym**4/2 + 3*s_sym*sigma_sym**6 - 11*s_sym*sigma_sym**4/4 + s_sym*sigma_sym**2,
            -7*s_sym**2*sigma_sym**6/8 + s_sym**2*sigma_sym**4/4 - 3*sigma_sym**8/4 + 3*sigma_sym**6/4,
            s_sym*sigma_sym**8/2 - s_sym*sigma_sym**6/8 - s_sym*sigma_sym**4/4,
            -s_sym**2*sigma_sym**6/8 - 3*sigma_sym**10/32 + sigma_sym**8/4 - sigma_sym**6/8
        ]
        D3 = sum(coeff * m1_sym**(8-i) for i, coeff in enumerate(D3_coeffs))
        
        # D4 expression (descending powers of m1)
        D4_coeffs = [
            -s_sym**2*sigma_sym**4/2 - 9*sigma_sym**4/4,
            3*s_sym**3*sigma_sym**4/2 + 27*s_sym*sigma_sym**4/4,
            -3*s_sym**4*sigma_sym**4/2 - 5*s_sym**2*sigma_sym**6/4 - 21*s_sym**2*sigma_sym**4/4 - 57*sigma_sym**6/8 + 27*sigma_sym**4/4,
            s_sym**5*sigma_sym**4/2 + 3*s_sym**3*sigma_sym**6 - 3*s_sym**3*sigma_sym**4/4 + 35*s_sym*sigma_sym**6/2 - 27*s_sym*sigma_sym**4/2,
            -9*s_sym**4*sigma_sym**6/4 + 3*s_sym**4*sigma_sym**4/2 - 5*s_sym**2*sigma_sym**8/4 - 45*s_sym**2*sigma_sym**6/4 + 21*s_sym**2*sigma_sym**4/4 - 141*sigma_sym**8/16 + 57*sigma_sym**6/4 - 27*sigma_sym**4/4,
            s_sym**5*sigma_sym**6/2 + 9*s_sym**3*sigma_sym**8/4 + 3*s_sym**3*sigma_sym**4/2 + 33*s_sym*sigma_sym**8/2 - 83*s_sym*sigma_sym**6/4 + 27*s_sym*sigma_sym**4/4,
            -9*s_sym**4*sigma_sym**8/8 + 7*s_sym**4*sigma_sym**6/8 - 5*s_sym**2*sigma_sym**10/8 - 61*s_sym**2*sigma_sym**8/8 + 11*s_sym**2*sigma_sym**6/2 + s_sym**2*sigma_sym**4/2 - 171*sigma_sym**10/32 + 159*sigma_sym**16/16 - 57*sigma_sym**6/8 + 9*sigma_sym**4/4,
            s_sym**5*sigma_sym**8/8 + 3*s_sym**3*sigma_sym**10/4 + 5*s_sym**3*sigma_sym**8/16 + s_sym**3*sigma_sym**6/4 + 27*s_sym*sigma_sym**10/4 - 17*s_sym*sigma_sym**8/2 + 13*s_sym*sigma_sym**6/4,
            -3*s_sym**4*sigma_sym**10/16 - 5*s_sym**2*sigma_sym**12/32 - 27*s_sym**2*sigma_sym**10/16 + 7*s_sym**2*sigma_sym**8/8 - s_sym**2*sigma_sym**6/8 - 51*sigma_sym**12/32 + 21*sigma_sym**10/8 - 9*sigma_sym**8/8,
            3*s_sym**3*sigma_sym**12/32 + s_sym**3*sigma_sym**10/16 - s_sym**3*sigma_sym**8/8 + 65*s_sym*sigma_sym**12/64 - 5*s_sym*sigma_sym**10/8 - 5*s_sym*sigma_sym**8/16,
            -s_sym**4*sigma_sym**10/32 - s_sym**2*sigma_sym**14/64 - s_sym**2*sigma_sym**12/64 - 5*s_sym**2*sigma_sym**10/32 - 3*sigma_sym**14/16 + 5*sigma_sym**12/16 - sigma_sym**10/8
        ]
        D4 = sum(coeff * m1_sym**(10-i) for i, coeff in enumerate(D4_coeffs))
        
        # Substitute m1 expression
        D1_sub = D1
        D2_sub = D2.subs(m1_sym, m1_expr)
        D3_sub = D3.subs(m1_sym, m1_expr)
        D4_sub = D4.subs(m1_sym, m1_expr)
        
        # Substitute s value (if provided)
        if s_val_str:
            D2_sub = D2_sub.subs(s_sym, s_val)
            D3_sub = D3_sub.subs(s_sym, s_val)
            D4_sub = D4_sub.subs(s_sym, s_val)
        
        # Simplify
        D2_sub = sp.simplify(D2_sub)
        D3_sub = sp.simplify(D3_sub)
        D4_sub = sp.simplify(D4_sub)
        
        return {
            'D1': D1_sub,
            'D2': D2_sub,
            'D3': D3_sub,
            'D4': D4_sub,
            'm1_expr': m1_expr,
            'success': True
        }
        
    except Exception as e:
        return {
            'error': str(e),
            'success': False
        }

def print_expressions_latex(expr_dict, m1_expr_str, s_val_str=None):
    """Print LaTeX formatted output"""
    if not expr_dict['success']:
        print(f"Error: {expr_dict['error']}")
        return
    
    print("=" * 80)
    if s_val_str:
        print(f"Input: m₁ = {m1_expr_str}, s = {s_val_str}")
    else:
        print(f"Input: m₁ = {m1_expr_str}")
    print("=" * 80)
    
    print("\nLaTeX output format:")
    print("-" * 40)
    
    # D1
    print(f"$$ D_1 = {sp.latex(expr_dict['D1'])} $$")
    
    # D2
    print(f"\n$$ D_2 = {sp.latex(expr_dict['D2'])} $$")
    
    # D3
    print(f"\n$$ D_3 = {sp.latex(expr_dict['D3'])} $$")
    
    # D4
    print(f"\n$$ D_4 = {sp.latex(expr_dict['D4'])} $$")
    
    print("\nSimplified output:")
    print("-" * 40)
    
    # Attempt to show polynomial expansion
    sigma_sym = sp.symbols('sigma')
    try:
        # D2 polynomial expansion
        if sigma_sym in expr_dict['D2'].free_symbols:
            poly2 = sp.Poly(expr_dict['D2'], sigma_sym)
            if poly2.degree() >= 0:
                print("\nD₂ polynomial expansion:")
                terms2 = []
                for i in range(poly2.degree(), -1, -1):
                    coeff = poly2.coeff_monomial(sigma_sym**i)
                    if coeff != 0:
                        terms2.append(f"{coeff}·σ^{i}")
                print("D₂ = " + " + ".join(terms2).replace("+ -", "- "))
    except:
        pass
    
    try:
        # D3 polynomial expansion
        if sigma_sym in expr_dict['D3'].free_symbols:
            poly3 = sp.Poly(expr_dict['D3'], sigma_sym)
            if poly3.degree() >= 0:
                print("\nD₃ polynomial expansion:")
                terms3 = []
                for i in range(poly3.degree(), -1, -1):
                    coeff = poly3.coeff_monomial(sigma_sym**i)
                    if coeff != 0:
                        terms3.append(f"{coeff}·σ^{i}")
                print("D₃ = " + " + ".join(terms3).replace("+ -", "- "))
    except:
        pass
    
    print("\n" + "=" * 80)

def interactive_calculator():
    """Interactive calculator"""
    print("D Expression Calculator")
    print("=" * 50)
    print("Available symbols: sigma (σ), s")
    print("Example inputs: sigma, sigma/2, -sigma^2, 1, 0.5*sigma")
    print("Type 'quit' to exit")
    print("=" * 50)
    
    while True:
        # Get user input
        m1_input = input("\nEnter m₁ expression: ").strip()
        
        if m1_input.lower() in ['quit', 'exit', 'q']:
            print("Goodbye!")
            break
        
        s_input = input("Enter s value (optional, press Enter to skip): ").strip()
        if not s_input:
            s_input = None
        
        # Compute and display
        expr_dict = get_D_expressions(m1_input, s_input)
        print_expressions_latex(expr_dict, m1_input, s_input)

# Directly run examples
if __name__ == "__main__":
    # Example 1: m₁ = σ/2, s = 1
    print("Example 1: m₁ = σ/2, s = 1")
    expr1 = get_D_expressions("sigma/2", "1")
    print_expressions_latex(expr1, "σ/2", "1")
    
    # Example 2: m₁ = σ, s = 1
    print("\n\nExample 2: m₁ = σ, s = 1")
    expr2 = get_D_expressions("sigma", "1")
    print_expressions_latex(expr2, "σ", "1")
    
    # Example 3: m₁ = -σ², s = 1
    print("\n\nExample 3: m₁ = -σ², s = 1")
    expr3 = get_D_expressions("-sigma**2", "1")
    print_expressions_latex(expr3, "-σ²", "1")
    
    # Start interactive mode
    print("\n" + "=" * 80)
    print("Starting interactive mode...")
    interactive_calculator()