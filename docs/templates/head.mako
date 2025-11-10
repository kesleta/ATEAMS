
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>

<script>
	$(document).ready(function() {
		const buttons = document.querySelectorAll(".collapsible");

		buttons.forEach(button => {
			var tableId = "#" + button.id + "-table";
			button.addEventListener("click", event => {
				$(tableId).toggle();
			});
		});
	})
</script>

<!-- Sidebar edits -->
<style>
	a:hover {
		text-decoration: underline;
	}

	a, a:visited {
		color: #058;
		transition: none;
	}

	header > h1 { display: none; }
	.toc { display: none; }

	img.resize {
		min-width: 10vw;
		max-width: 15vw;
		display: block;
		margin: 0 auto;
	}

	div.proof {
		border: 1px solid black;
		padding: 0em 1em;
		width: 90%;
		margin: 1em auto;
	}

	.tombstone {
		margin-top: -2em;
		float: right;
	}

	#index .two-column {
		column-count: 1 !important;
	}

	#ateams {
		margin-top: 0;
	}

	#sidebar {
		width: 20vw;
		font-size: 1rem;
		scrollbar-width: none;
	}
</style>

<!-- table edits -->
<style>
	.collapsible {
		border: 1px solid black;
		background: none;
		cursor: pointer;
		padding: 2.5px 5px;
		font-weight: bold;
	}

	.dataframe.profile-data-table {
		margin: auto;
		width: 75%;
		display: none;
	}

	.dataframe.profile-data-table thead tr {
		text-align: center !important;
	}

	.dataframe.profile-data-table thead {
		border-bottom: 3px solid black;
	}

	.dataframe.profile-data-table tbody tr th {
		vertical-align: middle;
	}

	.dataframe.profile-data-table tbody tr td {
		vertical-align: middle;
	}

	.dataframe.profile-data-table tbody > tr:nth-child(4n) {
		border-bottom: 1px solid black;
	}

	.dataframe.profile-data-table tbody td {
		text-align: center;
	}
</style>

<link rel="preload" as="image" href="https://mason.gmu.edu/~apizzime/assets/img/ic.4.cycle.jpg">
