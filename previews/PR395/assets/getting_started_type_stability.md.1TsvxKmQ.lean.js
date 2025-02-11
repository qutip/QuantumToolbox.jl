import{_ as t,c as e,o as p,ai as i,j as a,a as n}from"./chunks/framework.wlVbWh_x.js";const C=JSON.parse('{"title":"The Importance of Type-Stability","description":"","frontmatter":{},"headers":[],"relativePath":"getting_started/type_stability.md","filePath":"getting_started/type_stability.md","lastUpdated":null}'),l={name:"getting_started/type_stability.md"},h={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},o={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.332ex"},xmlns:"http://www.w3.org/2000/svg",width:"3.524ex",height:"2.732ex",role:"img",focusable:"false",viewBox:"0 -1060.7 1557.7 1207.4","aria-hidden":"true"},k={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},r={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.441ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 1079 1000","aria-hidden":"true"};function d(c,s,g,u,y,E){return p(),e("div",null,[s[18]||(s[18]=i("",37)),a("p",null,[s[2]||(s[2]=n("Before making a practical example, let's see the internal structure of the ")),s[3]||(s[3]=a("a",{href:"/QuantumToolbox.jl/previews/PR395/resources/api#QuantumToolbox.QuantumObject"},[a("code",null,"QuantumObject")],-1)),s[4]||(s[4]=n(" type. As an example, we consider the case of three qubits, and we study the internal structure of the ")),a("mjx-container",h,[(p(),e("svg",o,s[0]||(s[0]=[i("",1)]))),s[1]||(s[1]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("msubsup",null,[a("mrow",{"data-mjx-texclass":"ORD"},[a("mover",null,[a("mi",null,"σ"),a("mo",{stretchy:"false"},"^")])]),a("mi",null,"x"),a("mrow",{"data-mjx-texclass":"ORD"},[a("mo",{stretchy:"false"},"("),a("mn",null,"2"),a("mo",{stretchy:"false"},")")])])])],-1))]),s[5]||(s[5]=n(" operator:"))]),s[19]||(s[19]=i("",34)),a("p",null,[s[8]||(s[8]=n("In some functions of ")),s[9]||(s[9]=a("code",null,"QuantumToolbox.jl",-1)),s[10]||(s[10]=n(", you may find the use of the ")),s[11]||(s[11]=a("a",{href:"https://docs.julialang.org/en/v1/base/base/#Base.Val",target:"_blank",rel:"noreferrer"},[a("code",null,"Val")],-1)),s[12]||(s[12]=n(" type in the arguments. This is a trick to pass a value at compile time, and it is very useful to avoid type instabilities. Let's make a very simple example, where we want to create a Fock state ")),a("mjx-container",k,[(p(),e("svg",r,s[6]||(s[6]=[i("",1)]))),s[7]||(s[7]=a("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[a("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[a("mo",{"data-mjx-texclass":"ORD",stretchy:"false"},"|"),a("mi",null,"j"),a("mo",{fence:"false",stretchy:"false"},"⟩")])],-1))]),s[13]||(s[13]=n(" of a given dimension ")),s[14]||(s[14]=a("code",null,"N",-1)),s[15]||(s[15]=n(", and we give the possibility to create it as a sparse or dense vector. At first, we can write the function without using ")),s[16]||(s[16]=a("code",null,"Val",-1)),s[17]||(s[17]=n(":"))]),s[20]||(s[20]=i("",18))])}const F=t(l,[["render",d]]);export{C as __pageData,F as default};
