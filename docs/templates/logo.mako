

<!-- include a script for adding stuff to the end of proofs. -->
<script>
	// Get all the proofs in the document.
	proofs = document.getElementsByClassName("proof");

	// For each of the proofs, attach a floating child element in the bottom-right
	// corner.
	for (var proof of proofs) {
		// Create a proof-ending tombstone.
		square = document.createElement("div");
		square.className = "tombstone";
		square.innerHTML = "◼️";

		// Attach the tombstone to the proof.
		proof.appendChild(square);
	}
</script>


<header>
	<a class="homelink" rel="home" title="ATEAMS" href="https://github.com/apizzimenti/ATEAMS">
		<img class="resize" src="https://github.com/apizzimenti/ATEAMS/blob/feature/performant-sparse/docs/essential-cycle.jpeg?raw=true" alt="Homological percolation on the torus.">
	</a>
</header>