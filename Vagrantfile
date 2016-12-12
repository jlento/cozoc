Vagrant.configure("2") do |config|
    config.vm.box = "ubuntu/xenial64"
    config.vm.provider "virtualbox" do |v|
       v.memory = 3072
    end
    config.vm.provision "ansible_local" do |ansible|
        ansible.playbook = "playbook.yml"
    end
end
