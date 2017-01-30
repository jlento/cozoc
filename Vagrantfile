Vagrant.configure("2") do |config|
    config.vm.box = "ubuntu/xenial64"
    config.vm.provider "virtualbox" do |v|
      v.memory = 3072
      vb.customize ["guestproperty", "set", :id, "/VirtualBox/GuestAdd/VBoxService/--timesync-set-threshold", 10000]
    end
    config.vm.provision "ansible_local" do |ansible|
        ansible.playbook = "playbook.yml"
    end
end
