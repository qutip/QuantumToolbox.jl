import{_ as n,c as t,o as e,j as s,ai as l,a}from"./chunks/framework.DuQo_J4X.js";const o="/QuantumToolbox.jl/v0.25.2/assets/xghzttf.X7Reercr.svg",B=JSON.parse('{"title":"Solving for Steady-State Solutions","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/steadystate.md","filePath":"users_guide/steadystate.md","lastUpdated":null}'),h={name:"users_guide/steadystate.md"},r={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},p={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"6.599ex",height:"1.441ex",role:"img",focusable:"false",viewBox:"0 -626 2916.6 637","aria-hidden":"true"},d={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},k={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-1.577ex"},xmlns:"http://www.w3.org/2000/svg",width:"16.764ex",height:"4.939ex",role:"img",focusable:"false",viewBox:"0 -1486 7409.5 2183","aria-hidden":"true"},Q={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.489ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.618ex",height:"2.321ex",role:"img",focusable:"false",viewBox:"0 -810 1157.2 1026","aria-hidden":"true"},m={tabindex:"0"},T={style:{"text-align":"left"}},E={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},y={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.186ex"},xmlns:"http://www.w3.org/2000/svg",width:"6.979ex",height:"1.805ex",role:"img",focusable:"false",viewBox:"0 -716 3084.6 798","aria-hidden":"true"},u={style:{"text-align":"left"}},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},x={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.186ex"},xmlns:"http://www.w3.org/2000/svg",width:"6.979ex",height:"1.805ex",role:"img",focusable:"false",viewBox:"0 -716 3084.6 798","aria-hidden":"true"},v={style:{"text-align":"left"}},F={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},f={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.05ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.561ex",height:"1.645ex",role:"img",focusable:"false",viewBox:"0 -705 690 727","aria-hidden":"true"},b={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},w={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.05ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.561ex",height:"1.645ex",role:"img",focusable:"false",viewBox:"0 -705 690 727","aria-hidden":"true"},C={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},L={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.05ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.561ex",height:"1.645ex",role:"img",focusable:"false",viewBox:"0 -705 690 727","aria-hidden":"true"};function S(D,i,A,H,_,M){return e(),t("div",null,[i[43]||(i[43]=s("h1",{id:"doc:Solving-for-Steady-State-Solutions",tabindex:"-1"},[a("Solving for Steady-State Solutions "),s("a",{class:"header-anchor",href:"#doc:Solving-for-Steady-State-Solutions","aria-label":'Permalink to "Solving for Steady-State Solutions {#doc:Solving-for-Steady-State-Solutions}"'},"​")],-1)),i[44]||(i[44]=s("h2",{id:"introduction",tabindex:"-1"},[a("Introduction "),s("a",{class:"header-anchor",href:"#introduction","aria-label":'Permalink to "Introduction"'},"​")],-1)),s("p",null,[i[2]||(i[2]=a("For time-independent open quantum systems with decay rates larger than the corresponding excitation rates, the system will tend toward a steady state as ")),s("mjx-container",r,[(e(),t("svg",p,i[0]||(i[0]=[l("",1)]))),i[1]||(i[1]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"t"),s("mo",{stretchy:"false"},"→"),s("mi",{mathvariant:"normal"},"∞")])],-1))]),i[3]||(i[3]=a(" that satisfies the equation"))]),s("mjx-container",d,[(e(),t("svg",k,i[4]||(i[4]=[l("",1)]))),i[5]||(i[5]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("mfrac",null,[s("mrow",null,[s("mi",null,"d"),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"ρ"),s("mo",{stretchy:"false"},"^")])]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mtext",null,"ss")])])]),s("mrow",null,[s("mi",null,"d"),s("mi",null,"t")])]),s("mo",null,"="),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",{"data-mjx-variant":"-tex-calligraphic",mathvariant:"script"},"L")]),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"ρ"),s("mo",{stretchy:"false"},"^")])]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mtext",null,"ss")])]),s("mo",null,"="),s("mn",null,"0.")])],-1))]),s("p",null,[i[8]||(i[8]=a("Although the requirement for time-independence seems quite restrictive, one can often employ a transformation to the interaction picture that yields a time-independent Hamiltonian. For many these systems, solving for the asymptotic density matrix ")),s("mjx-container",Q,[(e(),t("svg",g,i[6]||(i[6]=[l("",1)]))),i[7]||(i[7]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"ρ"),s("mo",{stretchy:"false"},"^")])]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mtext",null,"ss")])])])],-1))]),i[9]||(i[9]=a(" can be achieved using direct or iterative solution methods faster than using master equation or Monte-Carlo simulations. Although the steady state equation has a simple mathematical form, the properties of the Liouvillian operator are such that the solutions to this equation are anything but straightforward to find."))]),i[45]||(i[45]=s("h2",{id:"Steady-State-solvers-in-QuantumToolbox.jl",tabindex:"-1"},[a("Steady State solvers in "),s("code",null,"QuantumToolbox.jl"),a(),s("a",{class:"header-anchor",href:"#Steady-State-solvers-in-QuantumToolbox.jl","aria-label":'Permalink to "Steady State solvers in `QuantumToolbox.jl` {#Steady-State-solvers-in-QuantumToolbox.jl}"'},"​")],-1)),i[46]||(i[46]=s("p",null,[a("In "),s("code",null,"QuantumToolbox.jl"),a(", the steady-state solution for a system Hamiltonian or Liouvillian is given by "),s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.steadystate"},[s("code",null,"steadystate")]),a(". This function implements a number of different methods for finding the steady state, each with their own pros and cons, where the method used can be chosen using the "),s("code",null,"solver"),a(" keyword argument.")],-1)),s("table",m,[i[31]||(i[31]=s("thead",null,[s("tr",null,[s("th",{style:{"text-align":"left"}},[s("strong",null,"Solver")]),s("th",{style:{"text-align":"left"}},[s("strong",null,"Description")])])],-1)),s("tbody",null,[s("tr",null,[i[17]||(i[17]=s("td",{style:{"text-align":"left"}},[s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.SteadyStateDirectSolver"},[s("code",null,"SteadyStateDirectSolver()")])],-1)),s("td",T,[i[12]||(i[12]=a("Directly solve ")),s("mjx-container",E,[(e(),t("svg",y,i[10]||(i[10]=[l("",1)]))),i[11]||(i[11]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"A"),s("mi",null,"x"),s("mo",null,"="),s("mi",null,"b")])],-1))]),i[13]||(i[13]=a(" using the standard method given in ")),i[14]||(i[14]=s("code",null,"Julia",-1)),i[15]||(i[15]=a()),i[16]||(i[16]=s("code",null,"LinearAlgebra",-1))])]),s("tr",null,[i[23]||(i[23]=s("td",{style:{"text-align":"left"}},[s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.SteadyStateLinearSolver"},[s("code",null,"SteadyStateLinearSolver()")])],-1)),s("td",u,[i[20]||(i[20]=a("Directly solve ")),s("mjx-container",c,[(e(),t("svg",x,i[18]||(i[18]=[l("",1)]))),i[19]||(i[19]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"A"),s("mi",null,"x"),s("mo",null,"="),s("mi",null,"b")])],-1))]),i[21]||(i[21]=a(" using the algorithms given in ")),i[22]||(i[22]=s("a",{href:"https://docs.sciml.ai/LinearSolve/stable/",target:"_blank",rel:"noreferrer"},[s("code",null,"LinearSolve.jl")],-1))])]),s("tr",null,[i[29]||(i[29]=s("td",{style:{"text-align":"left"}},[s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.SteadyStateEigenSolver"},[s("code",null,"SteadyStateEigenSolver()")])],-1)),s("td",v,[i[26]||(i[26]=a("Find the zero (or lowest) eigenvalue of ")),s("mjx-container",F,[(e(),t("svg",f,i[24]||(i[24]=[l("",1)]))),i[25]||(i[25]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",{"data-mjx-variant":"-tex-calligraphic",mathvariant:"script"},"L")])])],-1))]),i[27]||(i[27]=a(" using ")),i[28]||(i[28]=s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.eigsolve"},[s("code",null,"eigsolve")],-1))])]),i[30]||(i[30]=s("tr",null,[s("td",{style:{"text-align":"left"}},[s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.SteadyStateODESolver"},[s("code",null,"SteadyStateODESolver()")])]),s("td",{style:{"text-align":"left"}},[a("Solving time evolution with algorithms given in "),s("a",{href:"https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/",target:"_blank",rel:"noreferrer"},[s("code",null,"DifferentialEquations.jl"),a(" (ODE Solvers)")])])],-1))])]),i[47]||(i[47]=s("h2",{id:"Using-Steady-State-solvers",tabindex:"-1"},[a("Using Steady State solvers "),s("a",{class:"header-anchor",href:"#Using-Steady-State-solvers","aria-label":'Permalink to "Using Steady State solvers {#Using-Steady-State-solvers}"'},"​")],-1)),s("p",null,[i[36]||(i[36]=a("The function ")),i[37]||(i[37]=s("a",{href:"/QuantumToolbox.jl/v0.25.2/resources/api#QuantumToolbox.steadystate"},[s("code",null,"steadystate")],-1)),i[38]||(i[38]=a(" can take either a Hamiltonian and a list of collapse operators ")),i[39]||(i[39]=s("code",null,"c_ops",-1)),i[40]||(i[40]=a(" as input, generating internally the corresponding Liouvillian ")),s("mjx-container",b,[(e(),t("svg",w,i[32]||(i[32]=[l("",1)]))),i[33]||(i[33]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",{"data-mjx-variant":"-tex-calligraphic",mathvariant:"script"},"L")])])],-1))]),i[41]||(i[41]=a(" in Lindblad form, or alternatively, a Liouvillian ")),s("mjx-container",C,[(e(),t("svg",L,i[34]||(i[34]=[l("",1)]))),i[35]||(i[35]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",{"data-mjx-variant":"-tex-calligraphic",mathvariant:"script"},"L")])])],-1))]),i[42]||(i[42]=a(" passed by the user."))]),i[48]||(i[48]=l("",18))])}const V=n(h,[["render",S]]);export{B as __pageData,V as default};
