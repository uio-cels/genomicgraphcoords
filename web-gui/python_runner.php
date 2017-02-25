<?php
// Simple wrapper written in PHP for running python scripts through the webserver

header("Access-Control-Allow-Origin: *"); // Enable cors (so API can be called from anywhere)
error_reporting(E_ALL & ~E_NOTICE);

$request = $_SERVER['REQUEST_URI'];



$method = $_GET["method"];

if(isset($_GET["method"])){
	//$method = escapeshellarg($params[3]);
	if(!validate_string($method)) exit("Error #1");
	
	$output = array();
	$descriptorspec = array(
	   0 => array("pipe", "r"),  // stdin
	   1 => array("pipe", "w"),  // stdout
	   2 => array("pipe", "w"),  // stderr
	);

	//print_r($params);
	if(isset($_GET["params"])){
	        $p = str_replace(",", " ", $_GET["params"]);
	        $p = str_replace("%20", " ", $p);

			if(!validate_string($p)) exit("Error #2 .. $p");
			$arguments = $p; //str_replace(",", " ", $_GET["params"]);
	}
	else{
			$arguments = "";
	}
	$arguments = htmlspecialchars($arguments);

	//echo "<p>Arguments: " . $arguments . "</p>";
	//$command = "python3.4 /home/ivarandknut/checkout_genomic_intervals/interface.py $method $arguments";
	//$command = "python3 /home/ivarandknut/python-projects/OffsetBasedGraph/examples/gene_experiment.py $method $arguments";
	$command = "python3 /home/ivarandknut/checkout_genomic_intervals/gen_graph_coords.py $method $arguments";

	$process = proc_open($command, $descriptorspec, $pipes, dirname(__FILE__), null);

	$stdout = stream_get_contents($pipes[1]);
	fclose($pipes[1]);

	$stderr = stream_get_contents($pipes[2]);
	fclose($pipes[2]);

	if($method == ""){
			$stdout = explode("\n", $stdout);
	}
    //$stdout .= "command: $command";
	$out = array("stdout" => $stdout, "stderr" => $stderr);

	exit(json_encode($out));
}

echo "<p>Test sd</p>";
if(isset($_GET["make_dir"])){
	
	echo "Make dir"; 
	make_data_dir();
}

function make_data_dir(){
	$dir = "../data/tmp";
	try{
		echo "<p>Trying tmp</p>";
		$r = mkdir($dir, 0744);
		echo "<p>r: " . $r . "</p>";
		file_put_contents($dir.'/test.txt', 'Hello File');
	}
	catch(Exception $e){
		echo "<p>Error: " . $e->getMessage() ."</p>";
	}
	echo "MAde dir";
}

function validate_string($str) {
    return preg_match('/^[A-Za-z0-9_]+$/',$str);
}


?>