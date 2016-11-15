# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
  config.vm.box = "boxcutter/ubuntu1604-desktop"
  config.vm.provider "virtualbox" do |vb|
    vb.gui = true
    vb.cpus = 2
    vb.memory = 3072
  end
  config.vm.synced_folder ".", "/home/vagrant/cozoc", disabled: false
  config.vm.provision :shell, path: "bootstrap.sh"
end
