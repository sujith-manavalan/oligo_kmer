<!DOCTYPE html>
<html>
<head>
    <!-- add css style to the page -->
    <style>
        body {
            font-family: Helvetica;
        }

        h1 {
            font-size: 40px;
            color: #32c537;
	    margin-bottom: 8px;
	    text-align: center;       
 	}

        h3 {
            font-size: 20px;
            color: #942a2a;
            margin-bottom: 8px;
        }
    </style>
</head>
<body>
<h1>Oligo_kmer: K-mer based tool for the design of oligonucleotide probes</h1>

<!-- display details of the algorithm -->
<h3>Oligo-kmer algorithm highlights:</h3>
<ul>
  <li>Identify unique oligos sequences specific to each contig</li>
  <li>Hamming distance can be specified by the user to improve uniquness</li>
  <li>Oligo size can be specified by the user</li>
  <li>Remove oligos with repeat sequences</li>
  <li>Discard oligos with secondary structure</li>
  <li>Filter oligos based on GC content criteria</li>
  <li>Discard oligos with N</li>
  <li>Store the unique oligos in a SQL database</li>
</ul>

<!-- Process input from web user to collect the contig_id info -->
<form method="post" action="<?php echo $_SERVER['PHP_SELF'];?>">
  Enter the contig_id you are looking for: <input type="text" name="contig_id">
  <input type="submit">
</form>

<?php
/* variable to store contig_id */
$contig_id = "";

/* Checking the request method to see if the form was submitted */
if ($_SERVER["REQUEST_METHOD"] == "POST") {
    /* Collect user entered value */
    $contig_id = $_POST['contig_id'];
    echo "<br>";

    /* Check if user has entered contig_id, otherwise throw error */
    if (empty($contig_id)) {
        echo "No contig_id was provided please try again.<br>";
    } else {
        echo "<br>Contig_id provided was $contig_id.<br>";
    }
}

/* connect to the database using these details modify this parameters
to match the databse conbfiguration */
$server="localhost";
$username="smanavalan";
$password="";
$database="smanavalan";

/* establish a connection to the database */
$connect = mysqli_connect($server,$username,"",$database);

/* create an SQL query for fetching the unique kmers using 
FROM and WHERE statement */
$query = "SELECT kmer FROM unique_kmers WHERE contig_id = \"" . 
          $contig_id . "\" LIMIT 20";

/* Perform query using mysqli_query , store the result in $result */
$result = mysqli_query($connect, $query) 
          or trigger_error("Query Failed! SQL: $query - Error: "
          . mysqli_error($connect), E_USER_ERROR);

/* if statement checks if the query has any result , if blank if statement 
wont be executed */
if ($result = mysqli_query($connect, $query)) {
    /* fetch result from the row */
    echo "<br> Unique kmers(max 20 displayed) for 
          contig_id $contig_id:<br><br>";
    while ($row = mysqli_fetch_row($result)) {
        /* Print the unique kmer */
        printf(" %s<br>", $row[0]);
    }
    /* empty the memory using mysqli_free_result */
    mysqli_free_result($result);
} else {
    /* print an error message if no result */
    echo "No results";
}

/* close the database connection */
mysqli_close($connect);
?>
</body>
</html>
