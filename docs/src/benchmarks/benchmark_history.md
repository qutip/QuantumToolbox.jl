```@raw html
<iframe id="myIframe" src="https://albertomercurio.github.io/QuantumToolbox.jl/benchmarks/" style="width:100%; border:none; overflow:hidden;"></iframe>


<script>
  // Function to adjust the iframe height
  function adjustIframeHeight(event) {
    const iframe = document.getElementById('myIframe');
    if (event.origin === 'https://albertomercurio.github.io') { // Ensure this matches the iframe origin
      iframe.style.height = event.data + 'px';
    }
  }

  // Listen for messages from the iframe
  window.addEventListener('message', adjustIframeHeight, false);
</script>

```
