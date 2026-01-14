def get_samples_and_paths(sample_meth_beds, testing):     
    sample_ids = []
    meth_file_paths = []

    with open(sample_meth_beds, 'r') as file:
        for i, line in enumerate(file):
            if testing and i == 2: 
                print('[Testing] Using only 2 samples')
                break 

            # split() without arguments handles any whitespace (tabs or spaces)
            parts = line.strip().split()
            
            # Ensure the line isn't empty and has both parts
            if len(parts) == 2:
                sample_ids.append(parts[0])
                meth_file_paths.append(parts[1])

    if not testing: 
        print('[Production] Using all samples')

    print(f"Samples: {sample_ids}")

    return sample_ids, meth_file_paths

