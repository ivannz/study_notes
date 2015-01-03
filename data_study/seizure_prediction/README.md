<html>
<head>
	<title>American Epilepsy Society Seizure Prediction Challenge</title>
</head>
<body>
	<h2>Description</h2>
	<div>
		Intracranial EEG (iEEG) data clips are organized in folders containing training and testing data for each human or canine subject.<br/>
		The training data is organized into ten minute EEG clips labeled "Preictal" for pre-seizure data segments, or "Interictal" for non-seizure data segments.<br/>
		Each .mat file contains a data structure with fields as follow:
		<ol>
			<li>data: a matrix of EEG sample values arranged row x column as electrode x time.</li>
			<li>data_length_sec: the time duration of each data row</li>
			<li>sampling_frequency: the number of data samples representing 1 second of EEG data.</li>
			<li>channels: a list of electrode names corresponding to the rows in the data field</li>
			<li>sequence: the index of the data segment within the one hour series of clips. For example, preictal_segment_6.mat has a sequence number of 6, and represents the iEEG data from 50 to 60 minutes into the preictal data.</li>
		</ol>
		Preictal training and testing data segments are provided covering one hour prior to seizure with a five minute seizure horizon. (i.e. from 1:05 to 0:05 before seizure onset.) This pre-seizure horizon ensures that 1) seizures could be predicted with enough warning to allow administration of fast-acting medications, and 2) any seizure activity before the annotated onset that may have been missed by the epileptologist will not affect the outcome of the competition.<br/>
		Similarly, one hour sequences of interictal ten minute data segments are provided. The interictal data were chosen randomly from the full data record, with the restriction that interictal segments be as far from any seizure as can be practically achieved, to avoid contamination with preictal or postictal signals. In the long duration canine recordings it was possible to maintain a restriction of one week before or after a seizure. However, in the human recordings (which may be less than a week in total duration) interictal data was restricted to be more than four hours before or after any seizure.<br/>
	</div>
	<h2>Datasets</h2>
	<div>Datasets are stored on an external 1Tb storage hard drive in '/kaggle/eeg_data/'.</div>
</body>
</html>