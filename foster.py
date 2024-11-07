import streamlit as st
import schemdraw
import schemdraw.elements as elm
import sympy as sp

st.markdown("""
    <style>
        /* Set background and font color for the entire app */
        body {
            background-color: white;
            color: black;
        }
        .stApp {
            background-color: white;
            color: black;
        }
        
        /* Specific styling for titles, subtitles, and input fields */
        h1, h2, h3, h4, h5, h6, p, label, .stTextInput, .stSelectbox, .stMarkdown {
            color: black !important;
        }

        /* Styling for input text boxes */
        input, select {
            background-color: white !important;
            color: black !important;
        }

        /* Styling for button colors */
        .stButton>button {
            color: white;
            background-color: #0078d4;
        }

    </style>
""", unsafe_allow_html=True)

# Define the symbolic variable
s = sp.symbols('s')

def draw_Z_circuit(eqn, val):
    s = sp.symbols('s')
    print(val, s**val)
    mod_eqn=eqn/(s**val)
    foster_expansion = sp.apart(mod_eqn, s)

# Parse terms from the foster expansion
    terms = foster_expansion.as_ordered_terms()
    print("Terms", terms)

    for term in terms:
        term=term*(s**val)

    # with schemdraw.Drawing() as d:
    with schemdraw.Drawing() as d:
        d.config(bgcolor='white', color='black')
        for term in terms:
            mod_term=term*(s**val)
            print("Mod", mod_term)
            numerator, denominator = mod_term.as_numer_denom()
            print("Numerator", numerator, "denominator", denominator)
            if sp.degree(denominator) == 0 and numerator.has(s): # s implies an inductance in series
                L_val = term.coeff(s)
                d += elm.Inductor().label(f'L = {L_val} H')
            elif term.has(s**2):  # (a * s)/(d*(s**2 + b)) form implies parallel L and C
                # numerator, denominator = term.as_numer_denom()
                print(numerator, denominator)
                # Extract coefficients
                a = numerator.coeff(s)  # Coefficient of s in numerator
                d_val = denominator.coeff(s**2)  # Coefficient of s^2 in denominator
                b = denominator.subs(s, 0)/d_val
                print("a b d",a, d_val, b)
                L_val = a / (d_val * b)
                C_val = d_val / a
                print("L C", L_val, C_val)
                # Draw the parallel L and C components
                d.push()
                d += elm.Line().down().length(0.5)
                d += elm.Inductor().label(f'L = {L_val} H').down()
                # d += elm.Line().up().length(0.5)
                d.pop()
                d.push()
                d += elm.Line().right().length(3)
                d += elm.Line().down().length(0.5)
                d += elm.Capacitor().label(f'C = {C_val} F').down()
                d += elm.Line().left().length(3)
                d += elm.Line().right().length(3)
                # d.pop()
                # d += elm.Line().right()

            elif numerator.has(s) and sp.degree(denominator)==1: 
                # a*s/(s+b)
                a=numerator.coeff(s)/denominator.coeff(s)
                b=denominator.subs(s, 0)/denominator.coeff(s)
                L_val=a/b
                R_val=a
                d.push()
                d += elm.Line().down().length(0.5)
                d += elm.Inductor().label(f'L = {L_val} H').down()
                d += elm.Line().up().length(0.5)
                d.pop()
                d.push()
                d += elm.Line().right().length(3)
                d += elm.Line().down().length(0.5)
                d += elm.Resistor().label(f'R = {R_val} ohm').down()
                d += elm.Line().left().length(3)
                d += elm.Line().right().length(3)
                # d.pop()
                # d += elm.Line().right()

            elif denominator.has(s) and denominator.subs(s, 0)!= 0: # a/(s+b)
                a=numerator/denominator.coeff(s)
                b=denominator.subs(s, 0)/denominator.coeff(s)
                R_val=a/b
                C_val=1/a
                d.push()
                d += elm.Line().down().length(0.5)
                d += elm.Capacitor().label(f'C = {C_val} F').down()
                d += elm.Line().up().length(0.5)
                d.pop()
                d.push()
                d += elm.Line().right().length(3)
                d += elm.Line().down().length(0.5)
                d += elm.Resistor().label(f'R = {R_val} ohm').down()
                d += elm.Line().left().length(3)
                d += elm.Line().right().length(3)
                # d.pop()
            elif sp.degree(denominator) == 0 and sp.degree(numerator)==0:
                R_val = term
                d += elm.Resistor().label(f'R = {R_val} ohm')
    
            elif term.has(1/s):  # 1/s term implies a capacitance in series
                numerator, denominator = term.as_numer_denom()
                print("num den", numerator, denominator)
                # Extract coefficients
                m = denominator.coeff(s) 
                n = numerator

                C_val = m / n
                d += elm.Capacitor().label(f'C = {C_val} F')
            elif sp.degree(denominator) == 0 and sp.degree(numerator)==0:
                R_val = term
                d += elm.Resistor().label(f'R = {R_val} ohm')
            d += elm.Line().right()
        # Save the drawing to a file
        d.save('circuit.png')  # Save as 'circuit.png'
        
    # Display the image in Streamlit
    st.image('circuit.png')
        # d.show()

def draw_Y_circuit(eqn,val):
    s = sp.symbols('s')
    mod_eqn=eqn/(s**val)
    foster_expansion = sp.apart(mod_eqn, s)
    terms = foster_expansion.as_ordered_terms()
    with schemdraw.Drawing() as d:
        for term in terms:
            mod_term=term*(s**val)
            numerator, denominator = mod_term.as_numer_denom()
            if sp.degree(denominator) == 0 and numerator.has(s):  # s implies a capacitance in parallel
                C_val = 1 / term.coeff(s)
                d.push()
                d += elm.Line().right()
                d += elm.Line().up()
                d += elm.Capacitor().label(f'C = {C_val} F').up()
                d += elm.Line().left()
                d.pop()
            elif term.has(s**2):  # (a * s)/(d*(s**2 + b)) form implies parallel C and L
                # Extract coefficients
                a = numerator.coeff(s)  # Coefficient of s in numerator
                d_val = denominator.coeff(s**2)  # Coefficient of s^2 in denominator
                b = denominator.subs(s, 0)/d_val
                C_val = a / (b*d_val)
                L_val = d_val / a
                # Draw the parallel C and L components
                d.push()
                d += elm.Line().right()
                d += elm.Capacitor().label(f'C = {C_val} F').up()
                d += elm.Inductor().label(f'L = {L_val} H').up()
                d += elm.Line().left()
                d.pop()
                # d += elm.Line().right()

            elif numerator.has(s) and sp.degree(denominator) == 1: 
                # a*s/(s+b)
                a = numerator.coeff(s) / denominator.coeff(s)
                b = denominator.subs(s, 0) / denominator.coeff(s)
                C_val = a / b
                R_val = a
                d.push()
                d += elm.Line().right()
                d += elm.Capacitor().label(f'C = {C_val} F').up()
                d += elm.Resistor().label(f'R = {R_val} ohm').up()
                d += elm.Line().left()
                d.pop()
                # d += elm.Line().right()

            elif denominator.has(s) and denominator.subs(s, 0) != 0:  # a/(s+b)
                a = numerator / denominator.coeff(s)
                b = denominator.subs(s, 0) / denominator.coeff(s)
                R_val = a / b
                L_val = 1 / a
                d.push()
                d += elm.Line().right()
                d += elm.Resistor().label(f'R = {R_val} ohm').up()
                d += elm.Inductor().label(f'L = {L_val} H').up()
                d += elm.Line().left()
                d.pop()
                # d += elm.Line().right()

            # elif term.has(1/s):  # 1/s term implies an inductance in parallel
            #     numerator, denominator = term.as_numer_denom()
            #     m = denominator.coeff(s) 
            #     n = numerator

            #     L_val = n / m
            #     d.push()
            #     d += elm.Line().right()
            #     d += elm.Line().up()
            #     d += elm.Inductor().label(f'L = {L_val} H').up()
            #     d += elm.Line().left()
            #     d.pop()

            else:
                R_val = term
                d.push()
                d += elm.Line().right()
                d += elm.Line().up()
                d += elm.Resistor().label(f'R = {R_val} ohm').up()
                d += elm.Line().left()
                d.pop()
            d += elm.Line().right()

            
        # Save the drawing to a file
        d.save('circuit.png')  # Save as 'circuit.png'
        
    # Display the image in Streamlit
    st.image('circuit.png')

        # d.show() 
        #st.pyplot(d.fig)
        #st.plt(d.fig)

def find_typeofcircuit(function,foster_form):
    # get form input
    form=foster_form
    s = sp.symbols('s')
    numerator, denominator = function.as_numer_denom()
    zeros = sp.solve(numerator, s)
    poles = sp.solve(denominator, s)
    # print(poles)
    real_zeros = sorted([z for z in zeros if sp.im(z) == 0], key=lambda x: sp.re(x))
    real_poles = sorted([p for p in poles if sp.im(p) == 0], key=lambda x: sp.re(x))
   # Assuming zeros is already defined
    complex_zeros = sorted(
    [z for z in zeros if sp.im(z) != 0 or z == 0],  # Include 0 if it's present
    key=lambda x: (sp.re(x), sp.im(x))
)
    complex_poles = sorted(
    [p for p in poles if sp.im(p) != 0 or p == 0],  # Include 0 if it's present
    key=lambda x: (sp.re(x), sp.im(x))
)
    # print(complex_poles)
    alternating_real = []
    i, j = 0, 0  
    add_zero = real_zeros[0] < real_poles[0] if real_zeros and real_poles else True
    
    while i < len(real_zeros) or j < len(real_poles):
        if add_zero and i < len(real_zeros):
            alternating_real.append(real_zeros[i])
            i += 1
        elif not add_zero and j < len(real_poles):
            alternating_real.append(real_poles[j])
            j += 1
        add_zero = not add_zero

    alternating_complex = []
    i, j = 0, 0 
    add_zero = True if complex_zeros and not complex_poles else False
    if complex_zeros and complex_poles:
        add_zero = sp.re(complex_zeros[0]) < sp.re(complex_poles[0]) or \
                    (sp.re(complex_zeros[0]) == sp.re(complex_poles[0]) and sp.im(complex_zeros[0]) < sp.im(complex_poles[0]))

    while i < len(complex_zeros) or j < len(complex_poles):
        if add_zero and i < len(complex_zeros):
            alternating_complex.append(complex_zeros[i])
            i += 1
        elif not add_zero and j < len(complex_poles):
            alternating_complex.append(complex_poles[j])
            j += 1
        add_zero = not add_zero
    def is_sortedreal(arr):
        return all(arr[k] <= arr[k + 1] for k in range(len(arr) - 1))

    # Check if the alternating_real array is sorted
    is_alternating_real_sorted = is_sortedreal(alternating_real)

    valid_circuit=False

    if(is_alternating_real_sorted and (len(alternating_real) > 0) and alternating_real[-1]<=0):
        #print("its either RL or RC")
        if len(real_poles) > 0 and len(real_zeros) > 0:
         valid_circuit=True
         if(real_poles[-1]>real_zeros[-1]):
            print("this is an RC circuit")
            if(form==2):
                draw_Y_circuit(function,1)
            # return "RC"
         else:
            if(form==1):
                draw_Z_circuit(function,1)
            print("this is RL circuit")
            # return "RL"

    def is_sorted(arr):
        for k in range(len(arr) - 1):
            if sp.re(arr[k]) > sp.re(arr[k + 1]) or (sp.re(arr[k]) == sp.re(arr[k + 1]) and sp.im(arr[k]) > sp.im(arr[k + 1])):
                return False
        return True

    # Check if alternating_complex is sorted
    is_complex_sorted = is_sorted(alternating_complex)
    if is_complex_sorted and len(alternating_complex) != 0:
        valid_circuit=True
        print("This is an LC impedance circuit.")
        # return "LC"
    
    print(valid_circuit)
    if valid_circuit:
        print("here is a valid circuit")
        if form==1:
            draw_Z_circuit(function,0)
        elif form==2:
            draw_Y_circuit(function,0)

    # Display results
    print("Alternating Complex Array:", alternating_complex)
    print("Is the Alternating Complex Array Sorted?", is_complex_sorted)

    # Display results
    print("Real Zeros:", real_zeros)
    print("Real Poles:", real_poles)
    print("Complex Zeros:", complex_zeros)
    print("Complex Poles:", complex_poles)
    print("Alternating Real Array:", alternating_real)
    print("Alternating Complex Array:", alternating_complex)

# Main Streamlit app interface
st.title("Circuit Synthesis Web Application")

# User input for equation
# st.subheader("Input Transfer Function")
input_eqn = st.text_input("Enter transfer function (e.g., ((s + 1) * (s + 4))/(s*(s + 2) * (s + 5))):", "")
eqn = sp.sympify(input_eqn)

# Dropdown for choosing type of circuit
circuit_type = st.selectbox("Select circuit type:", ("Impedance (Z)", "Admittance (Y)"))
form=0
if(circuit_type=="Impedance (Z)"):
    form=1
else:
    form=2

find_typeofcircuit(eqn,form )


#---------------------------------------------------------


eqn=((s + 1) * (s + 4))/(s*(s + 2) * (s + 5))

# Z_s = ((s**2 + 2) * (s**2 + 4)) / (s * (s**2 + 3))
F_s = ((s**2 + 1) * (s**2 + 3)) / (s * (s**2 + 2)*(s**2 + 4))

#F_s=(s * (s**2 + 2))/((s**2 + 1) * (s**2 + 3))

#F_s=((s+1)*(s+3))/((s+2)*(s+6))

# Analyze the input and draw the circuit if equation is provided
# if input_eqn:
    # Parse the equation entered by the user
    # try:
    #     eqn = sp.sympify(input_eqn)
    #     st.write(f"Parsed Equation: {eqn}")

    #     # Button to start the drawing process
    #     if st.button("Draw Circuit Diagram"):
    #         if circuit_type == "Impedance (Z)":
    #             draw_circuit_diagram("Z", eqn, 1)
    #         else:
    #             draw_circuit_diagram("Y", eqn, 1)
    # except Exception as e:
    #     st.error(f"Error parsing equation: {e}")
