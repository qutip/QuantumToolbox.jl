import{_ as l,c as n,o as e,ai as t,j as s,a as i}from"./chunks/framework.DCbwEeOm.js";const _=JSON.parse('{"title":"Tensor Products and Partial Traces","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/tensor.md","filePath":"users_guide/tensor.md","lastUpdated":null}'),p={name:"users_guide/tensor.md"},h={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},o={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.357ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.395ex",height:"2.19ex",role:"img",focusable:"false",viewBox:"0 -810 1058.5 967.8","aria-hidden":"true"},r={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},d={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.357ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.224ex",height:"2.19ex",role:"img",focusable:"false",viewBox:"0 -810 982.8 967.8","aria-hidden":"true"},k={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},m={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.357ex"},xmlns:"http://www.w3.org/2000/svg",width:"7.555ex",height:"2.19ex",role:"img",focusable:"false",viewBox:"0 -810 3339.4 967.8","aria-hidden":"true"},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},Q={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"8.119ex",height:"1.971ex",role:"img",focusable:"false",viewBox:"0 -666 3588.6 871","aria-hidden":"true"},g={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},T={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.355ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.441ex",height:"1.358ex",role:"img",focusable:"false",viewBox:"0 -443 1079.1 600.1","aria-hidden":"true"},u={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},E={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.357ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.288ex",height:"1.359ex",role:"img",focusable:"false",viewBox:"0 -443 1011.2 600.8","aria-hidden":"true"},y={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},b={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.079ex",height:"1.464ex",role:"img",focusable:"false",viewBox:"0 -442 477 647","aria-hidden":"true"},C={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},x={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-1.552ex"},xmlns:"http://www.w3.org/2000/svg",width:"36.143ex",height:"4.084ex",role:"img",focusable:"false",viewBox:"0 -1119 15975.1 1805","aria-hidden":"true"};function F(v,a,w,f,L,B){return e(),n("div",null,[a[31]||(a[31]=t("",18)),s("p",null,[a[2]||(a[2]=i("To construct operators that act on an extended Hilbert space of a combined system, we similarly pass a list of operators for each component system to the ")),a[3]||(a[3]=s("a",{href:"/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.tensor"},[s("code",null,"tensor")],-1)),a[4]||(a[4]=i(" (or ")),a[5]||(a[5]=s("a",{href:"/QuantumToolbox.jl/v0.29.1/resources/api#Base.kron"},[s("code",null,"kron")],-1)),a[6]||(a[6]=i(") function. For example, to form the operator that represents the simultaneous action of the ")),s("mjx-container",h,[(e(),n("svg",o,a[0]||(a[0]=[t("",1)]))),a[1]||(a[1]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mi",null,"x")])])],-1))]),a[7]||(a[7]=i(" operator on two qubits:"))]),a[32]||(a[32]=t("",2)),s("p",null,[a[10]||(a[10]=i("To create operators in a combined Hilbert space that only act on a single component, we take the tensor product of the operator acting on the subspace of interest, with the identity operators corresponding to the components that are to be unchanged. For example, the operator that represents ")),s("mjx-container",r,[(e(),n("svg",d,a[8]||(a[8]=[t("",1)]))),a[9]||(a[9]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mi",null,"z")])])],-1))]),a[11]||(a[11]=i(" on the first qubit in a two-qubit system, while leaving the second qubit unaffected:"))]),a[33]||(a[33]=t("",5)),s("p",null,[a[16]||(a[16]=i("First, let’s consider a system of two coupled qubits. Assume that both the qubits have equal energy splitting, and that the qubits are coupled through a ")),s("mjx-container",k,[(e(),n("svg",m,a[12]||(a[12]=[t("",1)]))),a[13]||(a[13]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mi",null,"x")]),s("mo",null,"⊗"),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mi",null,"x")])])],-1))]),a[17]||(a[17]=i(" interaction with strength ")),s("mjx-container",c,[(e(),n("svg",Q,a[14]||(a[14]=[t("",1)]))),a[15]||(a[15]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"g"),s("mo",null,"="),s("mn",null,"0.05")])],-1))]),a[18]||(a[18]=i(" (in units where the bare qubit energy splitting is unity). The Hamiltonian describing this system is:"))]),a[34]||(a[34]=t("",7)),s("p",null,[a[25]||(a[25]=i("The simplest possible quantum mechanical description for light-matter interaction is encapsulated in the Jaynes-Cummings model, which describes the coupling between a two-level atom and a single-mode electromagnetic field (a cavity mode). Denoting the energy splitting of the atom and cavity ")),s("mjx-container",g,[(e(),n("svg",T,a[19]||(a[19]=[t("",1)]))),a[20]||(a[20]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"ω"),s("mi",null,"a")])])],-1))]),a[26]||(a[26]=i(" and ")),s("mjx-container",u,[(e(),n("svg",E,a[21]||(a[21]=[t("",1)]))),a[22]||(a[22]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"ω"),s("mi",null,"c")])])],-1))]),a[27]||(a[27]=i(", respectively, and the atom-cavity interaction strength ")),s("mjx-container",y,[(e(),n("svg",b,a[23]||(a[23]=[s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D454",d:"M311 43Q296 30 267 15T206 0Q143 0 105 45T66 160Q66 265 143 353T314 442Q361 442 401 394L404 398Q406 401 409 404T418 412T431 419T447 422Q461 422 470 413T480 394Q480 379 423 152T363 -80Q345 -134 286 -169T151 -205Q10 -205 10 -137Q10 -111 28 -91T74 -71Q89 -71 102 -80T116 -111Q116 -121 114 -130T107 -144T99 -154T92 -162L90 -164H91Q101 -167 151 -167Q189 -167 211 -155Q234 -144 254 -122T282 -75Q288 -56 298 -13Q311 35 311 43ZM384 328L380 339Q377 350 375 354T369 368T359 382T346 393T328 402T306 405Q262 405 221 352Q191 313 171 233T151 117Q151 38 213 38Q269 38 323 108L331 118L384 328Z",style:{"stroke-width":"3"}})])])],-1)]))),a[24]||(a[24]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"g")])],-1))]),a[28]||(a[28]=i(", the Jaynes-Cummings Hamiltonian can be constructed as:"))]),s("mjx-container",C,[(e(),n("svg",x,a[29]||(a[29]=[t("",1)]))),a[30]||(a[30]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("mi",null,"H"),s("mo",null,"="),s("mfrac",null,[s("msub",null,[s("mi",null,"ω"),s("mi",null,"a")]),s("mn",null,"2")]),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mi",null,"z")]),s("mo",null,"+"),s("msub",null,[s("mi",null,"ω"),s("mi",null,"c")]),s("msup",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"a"),s("mo",{stretchy:"false"},"^")])]),s("mo",null,"†")]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"a"),s("mo",{stretchy:"false"},"^")])]),s("mo",null,"+"),s("mi",null,"g"),s("mo",{stretchy:"false"},"("),s("msup",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"a"),s("mo",{stretchy:"false"},"^")])]),s("mo",null,"†")]),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mo",null,"−")]),s("mo",null,"+"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"a"),s("mo",{stretchy:"false"},"^")])]),s("msub",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mover",null,[s("mi",null,"σ"),s("mo",{stretchy:"false"},"^")])]),s("mo",null,"+")]),s("mo",{stretchy:"false"},")")])],-1))]),a[35]||(a[35]=t("",22))])}const D=l(p,[["render",F]]);export{_ as __pageData,D as default};
